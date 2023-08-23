import os
import glob
import re
import numpy as np
import pandas as pd
import string
import math
from scipy.stats import skew, skewtest

import cell_cycle_gating as ccg
from cell_cycle_gating import manual_gating as mg
from cell_cycle_gating import dead_cell_filter_ldrint as dcf_int

import patchworklib as pw
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.gridspec import GridSpec

import time
import multiprocessing

import plotnine
from plotnine import *
import seaborn as sns

from scipy.signal import find_peaks
from cell_cycle_gating import findpeaks as fp
from pomegranate.gmm import GeneralMixtureModel
from pomegranate.distributions import *
import torch

## Get the names of wells, valid wells (excluding outer two rows/columns), etc.
## input:
##    var: either "all_wells", "valid_wells", "all_well_rows", or "all_well_cols"
## output: a list with the corresponding info, e.g. ["A01", "A02", ... , "P24"] for all_wells
def get_well_names(var):
    all_well_rows = string.ascii_uppercase[0:16]
    all_well_cols = [str(num).zfill(2) for num in range(1,25)]
    valid_well_rows = all_well_rows[2:14]
    valid_well_cols = all_well_cols[2:22]
    ### get all wells
    all_wells = [row + col for row in all_well_rows for col in all_well_cols]
    all_wells.sort()
    ### get valid wells
    valid_wells = [row + col for row in valid_well_rows for col in valid_well_cols]
    valid_wells.sort()
    out = {
        'all_wells':all_wells,
        'valid_wells':valid_wells,
        'all_well_rows':all_well_rows,
        'all_well_cols':all_well_cols,
        'valid_well_rows':valid_well_rows,
        'valid_well_cols':valid_well_cols
    }
    return(out[var])

## Define a dictionary with the local folder name for experiments on each date
def define_folder_dict(name='folder_dict'):
    folder_dict = {
            '2020-11-17':'rep3',
            '2021-02-19':'rep4',
            '2021-02-26':'rep5',
            '2021-03-02':'rep6',
            '2021-04-06':'rep7',
            '2021-04-23':'rep8',
            '2021-05-18':'rep9',
            '2021-05-21':'rep10',
            '2021-06-11':'rep11',
            '2021-07-27':'2_rep1/210727-combo-rep1',
            '2021-07-30':'2_rep2/210730_combo_rep2',
            '2021-08-06':'2_rep3/210806_combo_rep3',
            '2021-10-05':'redo_rep1_and_2',
            '2021-10-15':'redo_rep1_and_2/redo_rep2',
            '2021-10-29':'redo_rep3'
        }
    globals()[name] = folder_dict
    return(None)

def get_base_dir():
    ### note: using linux/unix folder conventions -- would need to re-write for Windows
    base_dir = "/mnt/y/lsp-analysis/LINCS-combinations"
    if os.path.exists(base_dir):
        return(base_dir)
    elif os.path.exists("/Volumes/hits/lsp-analysis/LINCS-combinations"):
        return("/Volumes/hits/lsp-analysis/LINCS-combinations")
    else:
        raise Exception("Base directory not found -- need to mount research.files and supply its path, e.g. '/mnt/y/lsp-analysis/LINCS-combinations'")

def get_plates_to_regate():
    base_dir = get_base_dir()
    ### time_zero plates
    time_zero_file = os.path.join(base_dir, "re_gating", "plates_to_regate", "time_zero_regate_plates.csv")
    well_file = os.path.join(base_dir, "re_gating", "plates_to_regate", "time_zero_regate_wells.csv")
    df_plates1 = pd.read_csv(time_zero_file)
    df_wells1 = pd.read_csv(well_file)
    df_wells1['time'] = "time_zero"
    df_plates1['time'] = "time_zero"
    ### end-time control plates
    ctrl_end_file = os.path.join(base_dir, "re_gating", "plates_to_regate", "ctrl_end_time_regate_plates.csv")
    well_file = os.path.join(base_dir, "re_gating", "plates_to_regate", "ctrl_end_time_regate_wells.csv")
    df_plates2 = pd.read_csv(ctrl_end_file)
    df_wells2 = pd.read_csv(well_file)
    df_wells2['time'] = "end_time_control"
    df_plates2['time'] = "end_time_control"
    ### combine data frames
    df_plates = pd.concat([df_plates1, df_plates2])
    df_wells = pd.concat([df_wells1, df_wells2])
    return(df_plates, df_wells)

## Get a list of plate barcodes for a given date
## input:
##    date: e.g. '2021-10-15'
## output: list of barcodes e.g. '211015_combo_173'
def get_barcodes(date):
    date_formatted = date_format_switch(date)
    main_dir = get_data_dir(date = date)
    dirs = [ x for x in os.listdir(main_dir) if os.path.isdir( os.path.join(main_dir, x) )]
    ### match date at the start of the sub-directory
    dirs_barcodes = [ x for x in dirs if bool(re.match(date_formatted+"_combo", x)) ]
    return( dirs_barcodes )

## Switch the format of a date from YYYY-MM-DD to YYMMDD
## input:
##    date: e.g. '2021-02-19'
## output: e.g. '210219'
def date_format_switch(date):
    new_str = date[2:4] + date[5:7] + date[8:10]
    return(new_str)

## Get the date in YYYY-MM-DD format from a plate barcode
## input:
##    barcode: '210406_combo_71'
## output: e.g. '2021-04-06'
def date_from_barcode(barcode):
    date = '20' + barcode[0:2] + '-' + barcode[2:4] + '-' + barcode[4:6]
    return(date)

## Get the well-level data directory for a given date or barcode
## input:
##    barcode: a plate barcode, e.g. '210406_combo_71'
##    date: a date in YYYY-MM-DD format, e.g. '2021-04-06'
##    base_dir: the full path of the data folder, e.g. "/mnt/y/lsp-analysis/LINCS-combinations/"
## output:
##    returns the directory of the well-level data for a barcode
##    if a date is given and no barcode, returns the directory of all data for the date
##    if no date or barcode is given, returns the base directory of all data
def get_data_dir(barcode=None, date=None, base_dir = "/mnt/y/lsp-analysis/LINCS-combinations/"):
    ### note: using unix folder conventions -- would need to re-write for Windows
    ### set for osx
    if not os.path.exists(base_dir):
        base_dir = "/Volumes/hits/lsp-analysis/LINCS-combinations/"
    if barcode is None and date is None:
        return(base_dir)
    if date is None:
        date = date_from_barcode(barcode)
    if barcode is None:
        plate_dir = ''
    else:
        plate_dir = barcode
    #folder_dict = define_folder_dict()
    if not 'folder_dict' in globals(): define_folder_dict('folder_dict')
    local_dir = folder_dict[date]
    full_dir = os.path.join(base_dir, local_dir, plate_dir)
    return(full_dir)

## Get the filename for well-level intensities for a given barcode and well
## input:
##    barcode: a plate barcode, e.g. '210406_combo_71'
##    well: a well of interest, e.g. 'D06'
## output:
##    full path/filename of the well-level data
def get_well_file(barcode, well):
    date = date_from_barcode(barcode)
    data_dir = get_data_dir(barcode)
    ### example file style
    f1 = barcode+".result."+well+"[test].csv"
    #f2 = barcode+".result."+well+"[test].csv"
    files = os.listdir(data_dir)
    if f1 in files:
        well_file = os.path.join(data_dir, f1)
    else:
        print("well csv not found!")
    return(well_file)

## Read the well-level data for a given well and barcode
## input:
##    barcode: a plate barcode, e.g. '210406_combo_71'
##    well: a well of interest, e.g. 'D06'
## output:
##    a pandas dataframe of dye intensities for individual cells
def read_well_data(barcode, well):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    ff = get_well_file(barcode, well)
    df = pd.read_csv(ff)
    return(df)

## Re-name columns of well-level dataframe for LDR, DNA, EDU, etc.
## input: original data frame read from csv
## output: data frame with re-names columns
def rename_df_columns(df, silent=True, hoechst_as_dna=False):
    col_dict = {}
    ### check for well name
    if 'Well Name' in df.columns:
        if not silent: print("'Well Name' column found -- re-naming as 'well'")
        col_dict['Well Name'] = 'well'
    else:
        print(df.columns)
        if not silent: print('Well Name column not found')
        return(df)
    ### check for LDRint
    if 'ldrint' in df.columns:
        if not silent: print("'ldrint' column found -- re-naming as 'ldr'")
        col_dict['ldrint'] = 'ldr'
    else:
        print(df.columns)
        if not silent: print('ldrint column not found')
        return(df)
    ### check for DNAcontent/Hoechst
    dna_col1 = 'Cell: DNAcontent (DD-bckgrnd)'
    dna_col2 = 'Cell: DNAcontent (DDD-bckgrnd)'
    hoechst1 = 'Cell: HoechstINT (DD-bckgrnd)'
    hoechst2 = 'Cell: HoechstINT (DDD-bckgrnd)'
    check_dna1 = dna_col1 in df.columns
    check_dna2 = dna_col2 in df.columns
    check_hoechst1 = hoechst1 in df.columns
    check_hoechst2 = hoechst2 in df.columns
    if not (check_dna1 or check_dna2 or check_hoechst1 or check_hoechst2):
        if not silent: print(df.columns)
        if not silent: print('DNA column not found')
    else:
        ### if hoechst_as_dna, use the HoechstINT column as DNA if available
        if hoechst_as_dna:
            if check_hoechst1 and check_hoechst2:
                print('Warning: Two HoechstINT columns -- using ' + "'"+hoechst2+"'")
                dna_col = hoechst2
            elif check_hoechst1: dna_col = hoechst1
            elif check_hoechst2: dna_col = hoechst2
            elif check_dna1 and check_dna2:
                print('Warning: Two DNAcontent columns -- using ' + "'"+dna_col2+"'")
                dna_col = dna_col2
            elif check_dna1: dna_col = dna_col1
            elif check_dna2: dna_col = dna_col2
        ### otherwise use DNAConent column as DNA
        else:
            if check_dna1 and check_dna2:
                print('Warning: Two DNAcontent columns -- using ' + "'"+dna_col2+"'")
                dna_col = dna_col2
            elif check_dna1: dna_col = dna_col1
            elif check_dna2: dna_col = dna_col2
            elif check_hoechst1 and check_hoechst2:
                print('Warning: Two HoechstINT columns -- using ' + "'"+hoechst2+"'")
                dna_col = hoechst2
            elif check_hoechst1: dna_col = hoechst1
            elif check_hoechst2: dna_col = hoechst2
        if not silent: print("'"+dna_col+"'"+" column found -- re-naming as 'dna'")
        col_dict[dna_col] = 'dna'
    ### check for Edu (raw)
    if 'Cell: EdUrawINT (DDD-bckgrnd)' in df.columns:
        if not silent: print("'Cell: EdUrawINT (DDD-bckgrnd)' column found -- re-naming as 'edu_raw'")
        col_dict['Cell: EdUrawINT (DDD-bckgrnd)'] = 'edu_raw'
    elif 'Cell: EdUrawINT (DD-bckgrnd)' in df.columns:
        if not silent: print("'Cell: EdUrawINT (DD-bckgrnd)' column found -- re-naming as 'edu_raw'")
        col_dict['Cell: EdUrawINT (DD-bckgrnd)'] = 'edu_raw'
    ### check for Edu (background)
    if 'Cell: EdUbackground (DDD-bckgrnd)' in df.columns:
        if not silent: print("'Cell: EdUbackground (DDD-bckgrnd)' column found -- re-naming as 'edu_bg'")
        col_dict['Cell: EdUbackground (DDD-bckgrnd)'] = 'edu_bg'
    elif 'Cell: EdUbackground (DD-bckgrnd)' in df.columns:
        if not silent: print("'Cell: EdUbackground (DD-bckgrnd)' column found -- re-naming as 'edu_bg'")
        col_dict['Cell: EdUbackground (DD-bckgrnd)'] = 'edu_bg'
    ### re-name data-frame columns
    df = df.rename(columns=col_dict)
    if 'edu_raw' in df.columns and 'edu_bg' in df.columns:
        df['edu'] = df.edu_raw - df.edu_bg
    df['dna_colname'] = dna_col
    return(df)

def read_and_rename_well_data(barcode, well, silent = False, hoechst_as_dna=False):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    df = read_well_data(barcode, well)
    df = rename_df_columns(df, silent = silent, hoechst_as_dna=hoechst_as_dna)
    return(df)


## Get all wells with data for a given plate
## input:
##    barcode: a plate barcode, e.g. '210406_combo_71'
## output:
##    A list of wells, e.g. ['C03', 'C04', ... , 'N22']
def get_all_wells(barcode):
    data_dir = get_data_dir(barcode)
    wells_with_data = []
    all_wells = get_well_names("all_wells")
    for well in all_wells:
        ### example file style
        f1 = barcode+".result."+well+"[test].csv"
        f1_full = os.path.join(data_dir, f1)
        check = os.path.exists( f1_full )
        if check:
            wells_with_data.append(well)
    wells_with_data.sort()
    return(wells_with_data)

### maybe not necessary for LDR intensity data?
#def read_all_wells(barcode):
#    wells = get_all_wells(barcode)
#    df_list = [read_well_data(barcode, well) for well in wells]
#    return(df_list)

def get_ldr_cutoff(barcode, well, peak_loc = 1.2, silent=False, hoechst_as_dna=False):
    #df = read_well_data(barcode, well)
    #df = rename_df_columns(df)
    df = read_and_rename_well_data(barcode, well, silent, hoechst_as_dna=hoechst_as_dna)
    ldr_gates, ldr_lims = dcf_int.get_ldrgates(ldrint = df['ldr'], peak_loc=peak_loc) ## 1.2 is default
    return(ldr_gates[1])

def get_ldr_cutoff_many(barcode, wells, peak_loc = 1.2, silent=True, hoechst_as_dna=False):
    ldrs = [get_ldr_cutoff(barcode, well, peak_loc = peak_loc, silent=silent, hoechst_as_dna=hoechst_as_dna) for well in wells]
    return(ldrs)

def load_well_metadata(name = 'meta', folder=None, 
                       file='single_timepoint_cleaned_from_raw_2023-06-08.parquet'):
    ### read parquet file w/ all metadata
    if folder is None:
        folder = "/mnt/c/Users/NC168/git/LINCS_combos/data/cleaned/"
    full_file = os.path.join(folder, file)
    df = pd.read_parquet(full_file)
    globals()[name] = df
    return(None)

def get_wells(barcode, cell_line):
    ### get only the wells for a certain cell line on a given barcode
    if not 'meta' in globals(): load_well_metadata()
    #query = " cell_line == 'SUM1315' & barcode == '201117_combo_33' "
    query = "cell_line == '"+cell_line+"' & barcode == '"+barcode+"'"
    meta_sub = meta.query(query)
    wells = list(meta_sub.well)
    wells.sort()
    return(wells)

def get_cell_lines_on_plate(barcode):
    if not 'meta' in globals(): load_well_metadata()
    query = "barcode == '"+barcode+"'"
    meta_sub = meta.query(query)
    cell_lines = meta_sub.cell_line.unique()
    return(cell_lines)

def get_ldr_cutoffs_plate(barcode, peak_loc = 1.2, silent = True):
    cell_lines = get_cell_lines_on_plate(barcode)
    df = get_ldr_cutoffs_cell_line_and_barcode(barcode, cell_lines, peak_loc=peak_loc, silent=silent)
    return(df)

def get_ldr_cutoffs_cell_line_and_barcode(barcode, cell_lines, peak_loc=1.2, silent = True):
    df_list = []
    for cell_line in cell_lines:
        if not silent: print(cell_line)
        wells = get_wells(barcode, cell_line)
        ldrs = get_ldr_cutoff_many(barcode, wells, peak_loc=peak_loc)
        d = {'well':wells, 'ldr_cutoff': ldrs, 'barcode':barcode, 'cell_line': cell_line}
        df_tmp = pd.DataFrame(data=d)
        df_list.append(df_tmp)
    df = pd.concat(df_list)
    return(df)

## x_lims: tuple of x limits for the plot
## y_lims: tuple of y limits for the plot
def plot_ldr(df, peak_loc = 1.2, scatter = True, silent = True, show_fig = True, 
             fig=None, outer=None, i=None, title = "", x_lims=None, y_lims=None, add_ldr_line = None,
            ldr_gating_function = "new"):
    if ldr_gating_function == "new":
        res = get_ldrgates_new(ldrint = df['ldr'].copy())
        ldr_gates = res['ldr_gates']
    else:
        ldr_gates, _ = dcf_int.get_ldrgates(ldrint = df['ldr'].copy(), peak_loc=peak_loc)
    df = df.copy()
    ldr_cutoff = ldr_gates[1]
    #df = df.query("ldr > 0")
    #df['ldr'] = [x if x>0 else 10**(-10) for x in df.ldr]
    #### set negative ldr and dna values to the minimum positive values (just for plotting)
    df_pos1 = df.query("ldr > 0")
    min_ldr = np.min(df_pos1.ldr)
    df_pos2 = df.query("dna > 0")
    min_dna = np.min(df_pos2.dna)
    df['ldr'] = [x if x>0 else min_ldr for x in df.ldr]
    df['dna'] = [x if x>0 else min_dna for x in df.dna]
    if scatter:
        fig = mg.plot_ldr_dna_scatter(np.log10(df.dna), np.log10(df.ldr), ldr_cutoff, 
                                            dna_gates=None, plot_ldr_log10=True, is_ldrint=True,
                                           show_fig=show_fig, fig = fig, outer=outer, i=i,
                                     title=title, x_lims=x_lims, y_lims=y_lims, add_ldr_line = add_ldr_line)
    else:
        fig = mg.ldr_gating(np.log10(df.ldr), ldr_cutoff, nbins = 20)
    return(fig)

def plot_ldr_well(barcode, well, peak_loc = 1.2, scatter = True, silent = True, 
                  show_fig = True, fig=None, outer=None, i=None, title="", x_lims=None, 
                  y_lims=None, hoechst_as_dna=False, add_ldr_line=None, ldr_gating_function = "new"):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    df = read_and_rename_well_data(barcode, well, silent, hoechst_as_dna=hoechst_as_dna)
    #df = df.query("ldr > 0")
    #df['ldr'] = [x if x>0 else 10**(-10) for x in df.ldr]
    #df_pos1 = df.query("ldr > 0")
    #min_ldr = np.min(df_pos1.ldr)
    #df_pos2 = df.query("dna > 0")
    #min_dna = np.min(df_pos2.dna)
    #df['ldr'] = [x if x>0 else min_ldr for x in df.ldr]
    #df['dna'] = [x if x>0 else min_dna for x in df.dna]
    fig = plot_ldr(df, peak_loc = peak_loc, scatter = scatter, silent = silent, 
                   show_fig = show_fig, fig = fig, outer=outer, i=i, title = title,x_lims=x_lims,
                   y_lims=y_lims, add_ldr_line = add_ldr_line, ldr_gating_function = ldr_gating_function)
    return(fig)

def plot_ldr_many(barcode, wells, peak_loc = 1.2, scatter = False, silent = True, hoechst_as_dna=False):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    for well in wells:
        df = read_and_rename_well_data(barcode, well, silent, hoechst_as_dna=hoechst_as_dna)
        plot_ldr(df, peak_loc = peak_loc, scatter = scatter)

def plot_ldr_pdf(barcode, wells, peak_loc = 1.2, figname = "test_ldr.pdf", scatter = True, 
                 silent = True, show_fig = True, hoechst_as_dna=False):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    pdf_pages = PdfPages(figname)
    fig_list = []
    for i in range(len(wells)):
        well = wells[i]
        well_meta = get_well_meta(barcode, well)
        cell_line = list(well_meta.cell_line)[0]
        trt1 = list(well_meta.agent1)[0]
        trt2 = list(well_meta.agent2)[0]
        conc1 = list(well_meta.concentration1_chr)[0]
        conc2 = list(well_meta.concentration2_chr)[0]
        df = read_and_rename_well_data(barcode, well, silent, hoechst_as_dna=hoechst_as_dna)
        ldr_gates, ldr_lims = dcf_int.get_ldrgates(ldrint = df['ldr'], peak_loc=peak_loc)
        ldr_cutoff = ldr_gates[1]
        df = df.query("ldr > 0")
        #df['ldr'] = [x if x>0 else 10**(-10) for x in df.ldr]
        if scatter:
            fig = mg.plot_ldr_dna_scatter(np.log10(df.dna), np.log10(df.ldr), ldr_cutoff, dna_gates=None, 
                                          plot_ldr_log10=True, is_ldrint=True, show_fig = show_fig)
        else:
            fig = mg.ldr_gating(np.log10(df.ldr), ldr_cutoff, nbins = 20)
        fig_title = str(trt1) + ": "+ str(conc1) + " uM, " + str(trt2) + ": " + str(conc2) + " uM"
        #print(fig_title)
        fig.suptitle(well + "\n" + fig_title, fontsize=12)
        fig_list.append(fig)
        plt.close()
        pdf_pages.savefig(fig)
    pdf_pages.close()
    return(fig_list)

def test_regate(barcode, cell_line, peak_loc = 1.2, figname = "test_figure", scatter = True, silent = True, test = True, 
                show_fig = False, hoechst_as_dna=False):
    path1 = os.path.join('temp_regating', 'csv')
    path2 = os.path.join('temp_regating', 'pdf')
    if not os.path.exists(path1):
        os.makedirs(path1)
    if not os.path.exists(path2):
        os.makedirs(path2)
    wells = get_wells(barcode, cell_line)
    df_list = []
    csv_file = os.path.join(path1, figname+'.csv')
    pdf_file = os.path.join(path2, figname+'.pdf')
    plot_list = plot_ldr_pdf(barcode, wells, peak_loc, figname=pdf_file, scatter=scatter, silent=silent, 
                             show_fig=show_fig, hoechst_as_dna=hoechst_as_dna)
    print('figures written to: ' + pdf_file)
    for well in wells:
        df = read_and_rename_well_data(barcode, well, silent, hoechst_as_dna=hoechst_as_dna)
        df_tmp = dcf_int.get_counts_df(df=df, barcode=barcode, well=well, peak_loc = peak_loc)
        df_list.append(df_tmp)
    df_out = pd.concat(df_list)
    df_out.to_csv(csv_file)
    return(df_out, plot_list)

def plot_wells_ldr(barcode, cell_line, peak_loc=1.2, scatter = True, silent=True,
                   figname = None, output_dir="default_gating", hoechst_as_dna=False):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    if figname is None: figname=barcode+'_'+cell_line+'_'+'peak_loc_'+str(peak_loc)
    if not 'meta' in globals(): load_well_metadata()
    wells = get_wells(barcode, cell_line)
    df_list = []
    df_full_list = []
    for well in wells:
        df = read_and_rename_well_data(barcode, well, silent, hoechst_as_dna=hoechst_as_dna)
        #df = df.query("ldr > 0")
        df_tmp = dcf_int.get_counts_df(df=df, barcode=barcode, well=well, peak_loc = peak_loc)
        df_list.append(df_tmp)
        df_full_list.append(df)
    df2 = pd.concat(df_list)
    df_full = pd.concat(df_full_list)
    df_pos1 = df_full.query("ldr > 0")
    y_log = np.log10(df_pos1.ldr)
    y_lims = (min(y_log), max(y_log))
    print(y_lims)
    df_pos2 = df_full.query("dna > 0")
    x_log = np.log10(df_pos2.dna)
    x_lims = (min(x_log)-0.2, max(x_log)+0.2)
    ### save counts data frame to csv
    csv1 = "all_wells_" + figname + ".csv"
    csv1_full = os.path.join(output_dir, csv1)
    df2.to_csv(csv1_full)
    ### plot wells that changed
    fig_list = []
    pdf = "all_wells_scatter_" + figname + ".pdf"
    pdf_full = os.path.join(output_dir, pdf)
    pdf_pages = PdfPages(pdf_full)
    nb_plots = len(df2.well)
    plots_per_page = 6
    for i in range(nb_plots):
        #print(i)
        if i % plots_per_page == 0:
            fig = plt.figure(figsize=(8.5, 11))
            outer = GridSpec(3, 2, wspace=0.2, hspace=0.5)
        well = wells[i]
        #print(well)
        df = read_and_rename_well_data(barcode, well, silent=True, hoechst_as_dna=hoechst_as_dna)
        ### get well metadata
        well_meta = get_well_meta(barcode, well)
        cell_line = list(well_meta.cell_line)[0]
        trt1 = list(well_meta.agent1)[0]
        trt2 = list(well_meta.agent2)[0]
        conc1 = list(well_meta.concentration1_chr)[0]
        conc2 = list(well_meta.concentration2_chr)[0]
        ### add title to figures
        fig_title = str(trt1) + ": "+ str(conc1) + " uM, " + str(trt2) + ": " + str(conc2) + " uM"
        fig_title = well+", peak_loc = "+str(peak_loc)+"\n"+fig_title
    
        i_page = i % plots_per_page
        ### make figures
        fig_tmp = plot_ldr_well(barcode, well, peak_loc = peak_loc, scatter = scatter, 
                                 silent = silent, show_fig = False, fig = fig, outer = outer, i = i_page,
                               title = fig_title, x_lims=x_lims, y_lims=y_lims, hoechst_as_dna=hoechst_as_dna)
        #plt.close()
        fig_list.append(fig_tmp)
        if (i + 1) % plots_per_page == 0 or (i + 1) == nb_plots:
               plt.tight_layout()
               pdf_pages.savefig()
               plt.close('all')
    pdf_pages.close()
    return([df2, fig_list])

def get_well_meta(barcode, well):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    query = "barcode == '"+ barcode+ "' & well == '"+ well + "'"
    df_sub = meta.query(query)
    return(df_sub)

def plot_ldr_cutoff_change(barcode, cell_line, peak_loc, scatter = True, silent=True, 
                           default_peak_loc = 1.2, figname = None,
                          #output_dir="/mnt/y/lsp-analysis/LINCS-combinations/re_gating/new_gating"):
                           output_dir="temp_regating",
                          hoechst_as_dna=False):
    if figname is None: figname=barcode+'_'+cell_line+'_'+'peak_loc_'+str(peak_loc)
    if not 'meta' in globals(): load_well_metadata()
    wells = get_wells(barcode, cell_line)
    df_list_orig = []
    df_list_new = []
    for well in wells:
        df = read_and_rename_well_data(barcode, well, silent, hoechst_as_dna=hoechst_as_dna)
        #df = df.query("ldr > 0")
        df_tmp_new = dcf_int.get_counts_df(df=df, barcode=barcode, well=well, peak_loc = peak_loc)
        df_list_new.append(df_tmp_new)
        df_tmp_orig = dcf_int.get_counts_df(df=df, barcode=barcode, well=well, peak_loc = default_peak_loc)
        df_list_orig.append(df_tmp_orig)
    df_orig = pd.concat(df_list_orig)
    df_new = pd.concat(df_list_new)
    ### add suffixes to measured columns in each data frame
    df_orig2 = df_orig.rename(columns={c: c+'_orig' for c in df_orig.columns if c not in ['barcode', 'well']})
    df_new2 = df_new.rename(columns={c: c+'_new' for c in df_new.columns if c not in ['barcode', 'well']})
    ### join the data frames
    df2 = df_orig2.merge(df_new2, on = ['barcode', 'well'], how = 'inner')
    meta_select = meta[['barcode', 'cell_line', 'well', 'agent1', 'concentration1_chr', 'agent2', 'concentration2_chr', 'timepoint']]
    df2 = df2.merge(meta_select, on = ['barcode', 'well'], how = 'left')
    cols = ['barcode', 'well', 'cell_count__dead_orig', 'cell_count__dead_new', 'cell_count_orig', 'cell_count_new', 
                  'ldr_cutoff_orig', 'ldr_cutoff_new', 'cell_line','agent1', 'concentration1_chr', 'agent2', 
                  'concentration2_chr', 'timepoint']
    last_cols = [ x for x in df2.columns if x not in cols ]
    cols.extend(last_cols)
    df2 = df2[cols]
    df2['label'] = df2.apply(lambda row: row.well if row.ldr_cutoff_orig != row.ldr_cutoff_new else "", axis=1)
    df2_sub = df2.query("ldr_cutoff_orig != ldr_cutoff_new")
    df2_sub.reset_index(drop=True, inplace=True)
    ### plot new vs. old cutoffs, live/dead counts
    gg1 = pw.load_ggplot(ggplot(df2, aes(x = "ldr_cutoff_orig", y = "ldr_cutoff_new")) +\
        geom_point(alpha = 0.5) +\
        geom_label(aes(label="label"), alpha = 0.5, nudge_x = 0.05, nudge_y = 0.05), figsize=(3,3))
    gg2 = pw.load_ggplot(ggplot(df2, aes(x = "cell_count__dead_orig", y = "cell_count__dead_new")) +\
        geom_point(alpha = 0.5) +\
        geom_label(aes(label="label"), alpha = 0.5, nudge_x = 20, nudge_y = 5), figsize=(3,3))
    gg3 = pw.load_ggplot(ggplot(df2, aes(x = "cell_count_orig", y = "cell_count_new")) +\
        geom_point(alpha = 0.5) +\
        geom_label(aes(label="label"), alpha = 0.5, nudge_x = 100, nudge_y = 100), figsize=(3,3))
    
    gg = (gg1|gg2|gg3)
    pdf1 = figname + str("_summary.pdf")
    pdf1_full = os.path.join(output_dir, pdf1)
    gg.savefig(pdf1_full)
                               
    ### plot wells that changed
    fig_list_orig = []
    fig_list_new = []
    pdf2 = figname + str("_wells_changed.pdf")
    pdf2_full = os.path.join(output_dir, pdf2)
    pdf_pages = PdfPages(pdf2_full)
    nb_rows = len(df2_sub.well)
    rows_per_page = 3
    for i in range(nb_rows):
        if i % rows_per_page == 0:
            fig = plt.figure(figsize=(8.5, 11))
            outer = GridSpec(3, 2, wspace=0.2, hspace=0.5)
        well = df2_sub.well[i]
        df = read_and_rename_well_data(barcode, well, silent=True, hoechst_as_dna=hoechst_as_dna)
    
        ### get well metadata
        well_meta = get_well_meta(barcode, well)
        cell_line = list(well_meta.cell_line)[0]
        trt1 = list(well_meta.agent1)[0]
        trt2 = list(well_meta.agent2)[0]
        conc1 = list(well_meta.concentration1_chr)[0]
        conc2 = list(well_meta.concentration2_chr)[0]
        ### add title to figures
        fig_title = str(trt1) + ": "+ str(conc1) + " uM, " + str(trt2) + ": " + str(conc2) + " uM"
        fig_title_orig = well+", peak_loc = "+str(default_peak_loc)+" (default)"+"\n"+fig_title
        fig_title_new = well+", peak_loc = "+str(peak_loc)+"\n"+fig_title
    
        i_page = i % rows_per_page
        ### make figures for new and old peak_loc values
        fig_orig = plot_ldr_well(barcode, well, peak_loc = default_peak_loc, scatter = scatter, 
                                 silent = silent, show_fig = False, fig = fig, outer = outer, i = 2*i_page,
                                title = fig_title_orig, hoechst_as_dna=hoechst_as_dna)
        #plt.close()
        fig_new = plot_ldr_well(barcode, well, peak_loc = peak_loc, scatter = scatter, 
                                 silent = silent, show_fig = False, fig = fig, outer = outer, i = 2*i_page+1,
                               title = fig_title_new, hoechst_as_dna=hoechst_as_dna)
        #plt.close()
        fig_list_orig.append(fig_orig)
        fig_list_new.append(fig_new)
        if (i + 1) % rows_per_page == 0 or (i + 1) == nb_rows:
               plt.tight_layout()
               pdf_pages.savefig()
               plt.close('all')
    pdf_pages.close()
    ### write data frames to csv files
    # write cell counts for all wells
    csv1 = figname + str("_all_wells.csv")
    csv1_full = os.path.join(output_dir, csv1)
    df2.to_csv(csv1_full)
    # write cell counts for only wells where counts changed
    csv2 = figname + str("_wells_changed.csv")
    csv2_full = os.path.join(output_dir, csv2)
    df2_sub.to_csv(csv2_full)
    
    return([df2, df2_sub, gg, fig_list_orig, fig_list_new])

### plot LDR cutoffs
def plot_flagged_wells_ldr(barcode, cell_line, well_df, figname = None, peak_loc = 1.2, output_dir="default_gating", write_pdf=True,
                          hoechst_as_dna=False):
    if figname is None: figname=barcode+'_'+cell_line+'_'+'peak_loc_'+str(peak_loc)
    query = "cell_line == '"+cell_line+"' & barcode == '"+barcode+"'"
    well_df = well_df.query(query)
    wells = get_wells(barcode, cell_line)
    df_list = []
    for well in wells:
        df = read_and_rename_well_data(barcode, well, silent=True,hoechst_as_dna=hoechst_as_dna)
        #df = df.query("ldr > 0")
        df_tmp = dcf_int.get_counts_df(df=df, barcode=barcode, well=well, peak_loc = peak_loc)
        df_list.append(df_tmp)
    df = pd.concat(df_list)
    #return(df)
    df['flagged'] = ["flagged" if x in list(well_df.well) else "not_flagged" for x in df.well]
    gg = ggplot(df, aes(x = 'flagged', y = 'ldr_cutoff')) + geom_boxplot() + geom_jitter()

    csv_file = "flagged_wells_" + figname + ".csv"
    csv_full = os.path.join(output_dir, csv_file)
    pdf = "boxplot_ldr_cutoff_" + figname + ".pdf"
    pdf_full = os.path.join(output_dir, pdf)
    
    if write_pdf: gg.save(pdf_full, format = "pdf", width = 2.5, height = 3)
    well_df.to_csv(csv_full)
    return(df, gg)
    
def plot_problem_plate(barcode, cell_line, peak_loc=1.2, df_wells=None, scatter=True, silent=True, output_dir="default_gating",
                      hoechst_as_dna=False):
    final_dir = os.path.join(output_dir, cell_line + "_" + barcode)
    if not os.path.exists(final_dir): os.makedirs(final_dir)
    ### plot ldr vs. dna scatterplots for all wells:
    plot_wells_ldr(barcode, cell_line, peak_loc=peak_loc, scatter = scatter, silent=silent,
                   figname = None, output_dir=final_dir, hoechst_as_dna=hoechst_as_dna)
    ### plot ldr cutoffs for flagged vs. unflagged wells
    if df_wells is not None:
        plot_flagged_wells_ldr(barcode, cell_line, df_wells, output_dir = final_dir, write_pdf=True, hoechst_as_dna=hoechst_as_dna)

def plot_all_problem_plates(peak_loc=1.2, scatter=True, silent=True, output_dir = "default_gating", hoechst_as_dna=False):
    df_plates, df_wells = get_plates_to_regate()
    n_plates = len(df_plates.barcode)
    print("plotting LDR for " + str(n_plates) + " cell lines/plates")
    for i in range(n_plates):
        barcode = list(df_plates.barcode)[i]
        cell_line = list(df_plates.cell_line)[i]
        print("Plate " + str(i) + ": "+ cell_line + " " + barcode)
        plot_problem_plate(barcode, cell_line, peak_loc=peak_loc, df_wells=df_wells, scatter=scatter,
                           silent=silent, output_dir=output_dir, hoechst_as_dna=hoechst_as_dna)
        

### Todo:
## 1) def plot_problem_plate(): 
        ## steps:
        ## 1) call plot_wells_ldr to create plots of all wells for a problem plate-- write to pdf
        ## 2) load problem well metadata, create box plot of LDR cutoffs for "flagged wells" vs. "unflagged wells" -- print to pdf
        ## 3) save both (plus a csv, already written from plot_wells_ldr function) to "default_gating/
                               
## 2) loop over all bad plates/cell lines and plot LDR gating with default options
    ## for each, look at box plots and well-level scatter plots and pick a cutoff in-between the "flagged" and "non-flagged" ldr cutoffs.
    ## put new cutoffs manually into a dict/dataframe

## 3) call plot_LDR_cutoff_change with the new LDR cutoffs to write new plots and data to files

## Define a dictionary with the local folder name for experiments on each date
def define_regating_df(name = 'regate_df'):
    #barcode = '211015_combo_176'
    #cell_line = 'SUM1315'
    data = [
        ### Time zero plates
        ## note: plate 62 gating is fine with default peak_loc=1.2 on all but one well -- not sure why it was bad before
        {'barcode': '210406_combo_62', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False},
        {'barcode': '210423_combo_78', 'cell_line': 'SUM185PE', 'peak_loc': 2, 'hoechst_as_dna': False},
        ### End-time plates
        ### note: wells I06 and J07 -- almost 500 dead cells dead sub-g1, only ~50 LDR positive -- dna gating issue, not LDR gating?
        {'barcode': '210226_combo_51', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        ### note: only well I11 -- almost 500 dead cells dead sub-g1, only ~50 LDR positive -- dna gating issue, not LDR gating?
        ###  note: well I18 -- high subg1 dead cells
        {'barcode': '210226_combo_52', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        {'barcode': '210226_combo_53', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        {'barcode': '210226_combo_54', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        {'barcode': '210226_combo_55', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        {'barcode': '210226_combo_56', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        {'barcode': '210226_combo_57', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        {'barcode': '210302_combo_59', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        {'barcode': '210302_combo_60', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        {'barcode': '210302_combo_61', 'cell_line': 'HCC1937', 'peak_loc': 1.2, 'hoechst_as_dna': True}, #using Hoechst column as dna fixes dead count
        {'barcode': '210406_combo_69', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False}, #regating w/ default values gives low dead count
        {'barcode': '210406_combo_70', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False}, #regating w/ default values gives low dead count
        {'barcode': '210406_combo_71', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False}, #regating w/ default values gives low dead count
        # plate 72 -- E11 is the only control well that looks bad after regating -- will be solved by trimmed mean
        ### notes: a few non-control wells look wrong -- E05, possibly E07, F05, F19, etc.
        ### notes: a wide range of LDR cutoffs -- from 2 to 4 -- good cutoff looks like around 3 to 3.25
        {'barcode': '210406_combo_72', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False}, #regating w/ default values gives low dead count
        {'barcode': '210406_combo_73', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False}, #regating w/ default values gives low dead count
        {'barcode': '210406_combo_74', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False}, #regating w/ default values gives low dead count
        ### plate 75: E11 is the only control well that looks bad -- will be solved by trimmed mean
        {'barcode': '210406_combo_75', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False}, #regating w/ default values gives low dead count
        {'barcode': '210406_combo_76', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False}, #regating w/ default values gives low dead count
        {'barcode': '210406_combo_77', 'cell_line': 'SUM1315', 'peak_loc': 1.2, 'hoechst_as_dna': False}, #regating w/ default values gives low dead count
        {'barcode': '211005_combo_158', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211005_combo_160', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211005_combo_161', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211005_combo_162', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211005_combo_163', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211005_combo_164', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211005_combo_165', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211005_combo_166', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211015_combo_168', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211015_combo_169', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211015_combo_170', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211015_combo_171', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211015_combo_172', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211015_combo_173', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211015_combo_174', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211015_combo_175', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False},
        {'barcode': '211015_combo_176', 'cell_line': 'SUM1315', 'peak_loc': 2.75, 'hoechst_as_dna': False}
    ]
    df = pd.DataFrame(data)
    globals()[name] = df

def gate_well(barcode, well, peak_loc=1.2, silent=False, hoechst_as_dna=False):
    df = read_and_rename_well_data(barcode, well, silent, hoechst_as_dna=hoechst_as_dna)
    df_tmp = dcf_int.get_counts_df(df=df, barcode=barcode, well=well, peak_loc = peak_loc)
    return(df_tmp)

def regate_wells(silent=True):
    df_list2 = []
    for i in range(regate_df.shape[0]):
        print(i)
        barcode = regate_df['barcode'][i]
        cell_line = regate_df['cell_line'][i]
        peak_loc = regate_df['peak_loc'][i]
        hoechst_as_dna = regate_df['hoechst_as_dna'][i]
        if not 'meta' in globals(): load_well_metadata()
        wells = get_wells(barcode, cell_line)
        df_list = []
        for well in wells:
            #df = read_and_rename_well_data(barcode, well, silent, hoechst_as_dna=hoechst_as_dna)
            #df_tmp = dcf_int.get_counts_df(df=df, barcode=barcode, well=well, peak_loc = peak_loc)
            df_tmp = gate_well(barcode, well, peak_loc=peak_loc, silent=silent, hoechst_as_dna=hoechst_as_dna)
            df_list.append(df_tmp)
        df_out = pd.concat(df_list)
        df_list2.append(df_out)
    df_out2 = pd.concat(df_list2)
    return(df_out2)

def get_ldr_cutoffs_all(peak_loc = 1.2):
    if not 'folder_dict' in globals(): define_folder_dict('folder_dict')
    df_list_full = []
    for date in folder_dict.keys():
        print(date)
        plates = get_barcodes(date)
        df_list_date=[]
        for plate in plates:
            print(plate)
            df_tmp = get_ldr_cutoffs_plate(plate, peak_loc = peak_loc)
            df_list_date.append(df_tmp)
        df_date = pd.concat(df_list_date)
        df_list_full.append(df_date)
    df_full = pd.concat(df_list_full)
    return(df_full)

def get_ldr_cutoffs_fast(peak_loc=1.2):
    if not 'meta' in globals(): load_well_metadata()
    cutoffs = []
    for i in range(meta.shape[0]):
    #for i in range(50):
        if i % 1000 == 0: print(i)
        barcode = meta.barcode[i]
        well = meta.well[i]
        cutoff = get_ldr_cutoff(barcode, well, peak_loc = 1.2, silent=True)
        cutoffs.append(cutoff)
    return(cutoffs)

def get_ldr_cutoff_i(i, peak_loc=1.2):
    return(get_ldr_cutoff(meta.barcode[i], meta.well[i], peak_loc = 1.2, silent=True))

def get_ldr_cutoffs_parallel(peak_loc=1.2, nproc = 10, batch = 1000):
    if not 'meta' in globals(): load_well_metadata()
    cutoffs = []
    n_total = meta.shape[0]
    batches = np.ceil(meta.shape[0]/batch)
    for i in range(int(batches)):
        print(i)
        tic = time.time()
        start_batch = i*batch
        end_batch = min( (i+1)*batch, n_total)
        range_obj = range(start_batch, end_batch)
        pool = multiprocessing.Pool(nproc)
        cutoffs_batch = pool.map(get_ldr_cutoff_i, range_obj)
        #cutoffs_batch = pool.map(get_meta_i, range_obj)
        cutoffs.extend(cutoffs_batch)
        toc = time.time()
        print(str(toc-tic))
    return(cutoffs)

def get_counts_well(barcode, well, peak_loc=1.2, manual_ldr_cutoff=None, plot=True):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    df = read_and_rename_well_data(barcode, well)
    df_out = dcf_int.get_counts_df(df=df.copy(), barcode=barcode, well=well, 
                                   peak_loc = peak_loc, manual_ldr_cutoff=manual_ldr_cutoff)
    if plot: plot_ldr(df.copy(), peak_loc = peak_loc, add_ldr_line=manual_ldr_cutoff)
    return(df_out)

def kde_plot_wells(barcode, wells, title = "", add_legend=False, smoothing=1, column = "ldr"):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    fig, ax = plt.subplots()
    #ax.set_title(title, fontsize=12)
    for well in wells:
        df = read_and_rename_well_data(barcode, well, silent = True)
        ldrint = df[column].copy()
        ldrint[ldrint < 0] = float('nan')
        logint = np.log10(ldrint)
        #x, y = sns.kdeplot(logint, ax=ax, color = "grey", alpha=0.5).get_lines()[0].get_data()
        if add_legend:
            meta_sub=get_well_meta(barcode=barcode, well=well)
            trt = list(meta_sub.agent1)[0] + " " + list(meta_sub.concentration1_chr)[0] + "; " + list(meta_sub.agent2)[0] + " " + list(meta_sub.concentration2_chr)[0]
            sns.kdeplot(logint, ax=ax, label = well + " " + trt, alpha=0.25, bw_adjust=smoothing).set_title(title)
            ax.legend()
        else:
            sns.kdeplot(logint, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(title)
        plt.close()
    return(fig)

def get_kde_plot_data_well(barcode, well,smoothing=1):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    fig, ax = plt.subplots()
    df = read_and_rename_well_data(barcode, well, silent = True)
    ldrint = df['ldr'].copy()
    ldrint[ldrint < 0] = float('nan')
    logint = np.log10(ldrint)
    x, y = sns.kdeplot(logint, ax=ax, color = "grey", alpha=0.5, bw_adjust=smoothing).get_lines()[0].get_data()
    plt.close()
    return(x,y)

def get_logldrint(barcode, well, remove_na=False, return_ldrint_also=False):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    df = read_and_rename_well_data(barcode, well, silent = True)
    ldr_int_orig = df['ldr'].copy()
    ldrint = df['ldr'].copy()
    if not remove_na:
        ldrint[ldrint < 0] = float('nan')
        logint = np.log10(ldrint)
    else:
        ldrint = ldrint[ldrint > 0]
        logint = np.log10(ldrint)
        logint = logint[ [not x for x in np.isnan(logint)] ]
        logint = logint[ [not x for x in np.isinf(logint)] ]
    if return_ldrint_also:
        return(list(logint), ldr_int_orig)
    else:
        return(list(logint))

def kde_plot_cell_line(barcode, cell_line, n_wells=None, well_start=0, title = "", 
                       add_legend=False, smoothing=1, column="ldr"):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    wells = get_wells(barcode, cell_line)
    print(wells)
    if n_wells is not None: wells = wells[well_start:(well_start+n_wells)]
    if title == "": title = cell_line + " " + barcode
    fig = kde_plot_wells(barcode, wells, title = title, add_legend=add_legend, 
                         smoothing=smoothing, column = column)
    return(fig)

def kde_plot_avg(barcode, cell_line, smoothing=1):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    wells = get_wells(barcode, cell_line)
    #if n_wells is not None: wells = wells[0:n_wells]
    logint_all = []
    #ax = axs[row,col]
    fig, ax = plt.subplots()
    for well in wells:
        logint = get_logldrint(barcode, well)
        logint_all.extend(logint)
    #return(logint_all)
    sns.kdeplot(logint_all, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(barcode + " " + cell_line)
    plt.close()
    return(fig)

def kde_plot_all_avg(cell_line, n_barcodes=None, smoothing=1):
    barcodes = get_all_plates_for_cell_line(cell_line)
    fig, ax = plt.subplots()
    if n_barcodes is not None: barcodes = barcodes[0:n_barcodes]
    print(len(barcodes))
    for i in range(len(barcodes)):
        #print(i)
        barcode = barcodes[i]
        wells = get_wells(barcode, cell_line)
        #if n_wells is not None: wells = wells[0:n_wells]
        logint_all = []
        for well in wells:
            logint = get_logldrint(barcode, well)
            logint_all.extend(logint)
        #return(logint_all)
        sns.kdeplot(logint_all, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(cell_line)
        plt.close()
    return(fig)

def kde_plot_plate(barcode, n_wells = None, output_dir="", filename="test_kde.pdf",
                   add_ldr_line=False, add_median_ldr_line=False, smoothing=1):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    cell_lines = get_cell_lines_on_plate(barcode)
    pdf_full = os.path.join(output_dir, filename)
    nb_plots = len(cell_lines)
    if nb_plots != 6: print(barcode + ": " + str(nb_plots) + " cell lines")
    if nb_plots == 0: return(None)
    ncols = 2
    nrows = int(np.ceil(nb_plots/2))
    nrows = 3
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(9, 4*nrows),
                        layout="constrained", sharex= "all")
    ldr_cutoffs = []
    for row in range(nrows):
        for col in range(ncols):
            i = row*ncols + col
            if i == nb_plots: break
            cell_line = cell_lines[i]
            wells = get_wells(barcode, cell_line)
            if n_wells is not None: wells = wells[0:n_wells]
            logint_all = []
            ax = axs[row,col]
            for well in wells:
                logint = get_logldrint(barcode, well)
                logint_all.extend(logint)
                sns.kdeplot(logint, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(cell_line)
                ldr_gates, ldr_lims = dcf_int.get_ldrgates(np.array([10**x for x in logint]))
                ldr_cutoff = ldr_gates[1]
                ldr_cutoffs.append(ldr_cutoff)
                if add_ldr_line:
                    ax.axvline(ldr_cutoff,ymin=0, ymax=0.1, color = "red", linestyle = "--")
            #sns.kdeplot(logint_all, ax=ax, color = "red").set_title(cell_line)
            xmin = -2
            xmax = 6
            x_lims = (xmin, xmax)
            plt.xlim(x_lims)
            x_ticks = np.arange(np.ceil(xmin), np.floor(xmax)+1)
            plt.xticks(x_ticks)
            ax.tick_params(labelbottom=True)
            if add_median_ldr_line:
                med_cutoff = np.median(ldr_cutoffs)
                #print(med_cutoff)
                ax.axvline(x=med_cutoff, ymin=0, ymax=1, color = "orange")
    fig.suptitle(barcode)
    plt.savefig(pdf_full)
    return(None)

def kde_plot_all_plates(smoothing=1):
    if not 'folder_dict' in globals(): define_folder_dict('folder_dict')
    for date in folder_dict.keys():
        print(date)
        plates = get_barcodes(date)
        plates.sort()
        for plate in plates:
            folder = "density_plots"
            if not os.path.exists(folder): os.makedirs(folder)
            file = plate + ".pdf"
            if os.path.exists(os.path.join(folder, file)):
                print(plate)
                print("pdf already written")
            else:
                kde_plot_plate(plate, output_dir = folder, filename = file, smoothing=smoothing)
                print(plate)
    return(None)

def kde_plot_all_plates_median(smoothing=1):
    if not 'folder_dict' in globals(): define_folder_dict('folder_dict')
    for date in folder_dict.keys():
        print(date)
        plates = get_barcodes(date)
        plates.sort()
        for plate in plates:
            folder = "density_plots_median_ldr"
            if not os.path.exists(folder): os.makedirs(folder)
            file = plate + "_median.pdf"
            if os.path.exists(os.path.join(folder, file)):
                print(plate)
                print("pdf already written")
            else:
                kde_plot_plate(plate, output_dir = folder, filename = file, add_median_ldr_line=True, smoothing=smoothing)
                print(plate)
    return(None)

def get_all_plates_for_cell_line(cell_line):
    if not 'meta' in globals(): load_well_metadata('meta')
    query = "cell_line == '"+cell_line+"'"
    meta_sub = meta.query(query)
    barcodes = meta_sub.barcode.unique()
    return(barcodes)

def kde_plot_all_plates_cell_line(cell_line, n_wells = None, n_barcodes=None, output_dir="", 
                                  filename="test_kde_cell_line.pdf", smoothing=1):
    barcodes = get_all_plates_for_cell_line(cell_line)
    barcodes.sort()
    pdf_full = os.path.join(output_dir, filename)
    ncols = 2
    nrows = 3
    pdf_pages = PdfPages(pdf_full)
    fig_list = []
    nb_plots = len(barcodes)
    if n_barcodes is not None: nb_plots = n_barcodes
    for i in range(nb_plots):
        barcode = barcodes[i]
        i_rel = (i % (ncols*nrows))
        col = i_rel % ncols
        row = i_rel // ncols
        print(str(row) + ", " + str(col))
        if i % (ncols*nrows) == 0:
            fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(9, 12),
                        layout="constrained", sharex= "all")
            fig.suptitle(cell_line)
        wells = get_wells(barcode, cell_line)
        if n_wells is not None: wells = wells[0:n_wells]
        ax = axs[row,col]
        for well in wells:
            logint = get_logldrint(barcode, well)
            #logint_all.append(logint)
            sns.kdeplot(logint, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(barcode)
        #sns.kdeplot(logint_all, ax=ax, color = "red").set_title(cell_line)
        xmin = -2
        xmax = 6
        x_lims = (xmin, xmax)
        plt.xlim(x_lims)
        x_ticks = np.arange(np.ceil(xmin), np.floor(xmax)+1)
        plt.xticks(x_ticks)
        ax.tick_params(labelbottom=True)
        plt.close()
        if (i + 1) % (ncols*nrows) == 0 or (i + 1) == nb_plots:
            plt.tight_layout()
            pdf_pages.savefig(fig)
            plt.close('all')
    pdf_pages.close()
    return(None)

def get_mixture_cutoff(model, n_iter=20, prob_tol=0.0001, update_tol=0.0001, silent=True):
    mean1 = model.distributions[0].means.item()
    mean2 = model.distributions[1].means.item()
    #print(mean1)
    #print(mean2)
    min = np.min([mean1, mean2])
    max = np.max([mean1, mean2])
    prob = 0
    cutoff = (min+max)/2
    error = 0.5-prob
    i=0
    diff = 1
    while (abs(error) > prob_tol or diff > update_tol) and i<n_iter:
        i=i+1
        if not silent:
            print(i)
            print(cutoff)
        if i>1:
            last_cutoff=cutoff
        else:
            last_cutoff=min
        prob_both=model.predict_proba([[cutoff]])
        prob1 = prob_both[0][0].item()
        prob2 = prob_both[0][1].item()
        if mean1<mean2:
            prob=prob1
        else:
            prob=prob2
        error = 0.5-prob
        if prob>0.5:
            min=cutoff
        else:
            max=cutoff
        cutoff=(min+max)/2
        diff=last_cutoff-cutoff
        if not silent:
            print(prob_both)
            print(prob)
    return(cutoff)

### inputs: x,y from sns.kdeplot output
def get_peaks_ldr(x, y, smoothing=1, first_peak_min=0.5,
                 min_prominence=0, min_peak_height=0.02, min_peak_distance=0.5,
                 single_peak_cutoff=3, silent=True):
    x=x.copy()
    y=y.copy()
    peak_locs, peak_props = find_peaks(y.copy(), height=0, prominence=0)
    peak_ldrs = x[peak_locs]
    peak_props_orig=peak_props.copy()
    peak_ldrs_orig=peak_ldrs.copy()
    ### get rid of invalid peaks
    #print(peak_ldrs)
    #print(peak_props)
    keep_peaks = [peak_ldrs[i] > first_peak_min and 
                  peak_props['prominences'][i]>min_prominence and
                 peak_props['peak_heights'][i]>min_peak_height for i in range(len(peak_ldrs))]
    valid_peaks = np.sum(keep_peaks)
    indexes = [i for i in range(len(keep_peaks)) if keep_peaks[i]]
    peak_locs = [peak_locs[ind] for ind in indexes]
    peak_ldrs = [peak_ldrs[ind] for ind in indexes]
    for prop in peak_props:
        peak_props[prop] = np.array([ peak_props[prop][ind] for ind in indexes ])
    ### identify main peak and secondary peak
    peak_order = np.argsort(-np.array(peak_props['peak_heights'])) ### in order of peak height
    main_peak_index=peak_order[0]
    secondary_peak_index = None
    #### find the secondary peak
    if len(peak_locs)>1:
        ### loop over rest of peaks in order of height
        for i in range(1,len(peak_locs)):
            ### select next highest peak as long as it is sufficiently far from the first peak
            #if peak_ldrs[main_peak_index] < single_peak_cutoff:
            #    check1 = peak_ldrs[peak_order[i]] - peak_ldrs[main_peak_index] > min_peak_distance
            #else:
            #    check1 = peak_ldrs[main_peak_index] - peak_ldrs[peak_order[i]] > min_peak_distance
            check1 = abs(peak_ldrs[peak_order[i]] - peak_ldrs[main_peak_index]) > min_peak_distance
            if check1:
                secondary_peak_index = peak_order[i]
                break
                
    if secondary_peak_index is None:
        if not silent: print("search for shelf")
        ### search for a shelf if no secondary peak found
        shelf_dict = find_shelf(x.copy(), y.copy(), main_peak_ldr=peak_ldrs[main_peak_index], 
                                min_peak_height=min_peak_height,first_peak_min=first_peak_min)
        indexes = [main_peak_index]
        peak_locs = [peak_locs[ind] for ind in indexes]
        peak_ldrs = [peak_ldrs[ind] for ind in indexes]
        for prop in peak_props:
            peak_props[prop] = np.array([ peak_props[prop][ind] for ind in indexes ])
        ### add shelf if found
        if shelf_dict is not None:
            shelf_ldr = shelf_dict['ldr']
            shelf_height = shelf_dict['height']
            if shelf_ldr < peak_ldrs[0]:
                peak_ldrs = [shelf_ldr, peak_ldrs[0]]
                shelf = "peak1"
                peak_height1=shelf_height
                peak_height2=peak_props['peak_heights'][0]
            else:
                peak_ldrs = [peak_ldrs[0], shelf_ldr]
                shelf = "peak2"
                peak_height1=peak_props['peak_heights'][0]
                peak_height2=shelf_height
        else: ### single peak case
            shelf="none"
            if peak_ldrs[0]<single_peak_cutoff:
                peak_height1=peak_props['peak_heights'][0]
                peak_height2=float('nan')
            else:
                peak_height1=float('nan')
                peak_height2=peak_props['peak_heights'][0]
    else:
        shelf="none"
        ### add peaks in order of their ldr intensity
        if peak_ldrs[main_peak_index] < peak_ldrs[secondary_peak_index]:
            indexes = [main_peak_index, secondary_peak_index]
        else:
            indexes = [secondary_peak_index, main_peak_index]
        peak_locs = [peak_locs[ind] for ind in indexes]
        peak_ldrs = [peak_ldrs[ind] for ind in indexes]
        for prop in peak_props:
            peak_props[prop] = np.array([ peak_props[prop][ind] for ind in indexes ])
        peak_height1=peak_props['peak_heights'][0]
        peak_height2=peak_props['peak_heights'][1]
    if(len(peak_ldrs)==0):
        print("Warning! No valid peaks")
        shelf="none"
        out={'peak1': float('nan'), 'peak2': float('nan'), 'peak1_height':float('nan'), 'peak2_height':float('nan'),
             'shelf':shelf, 'peak_props':peak_props, 'peak_props_orig':peak_props_orig, 'peak_ldrs_orig': peak_ldrs_orig}
    if(len(peak_ldrs)==1):
        if(peak_ldrs[0] < single_peak_cutoff):
            if not silent: print("only live cell peak found")
            out={'peak1': peak_ldrs[0], 'peak2': float('nan'), 'peak1_height':peak_height1, 'peak2_height':peak_height2,
                 'shelf':shelf, 'peak_props':peak_props, 'peak_props_orig':peak_props_orig, 'peak_ldrs_orig': peak_ldrs_orig}
        else:
            if not silent: print("only dead cell peak found")
            out={'peak1': float('nan'), 'peak2': peak_ldrs[0],  'peak1_height':peak_height1, 'peak2_height':peak_height2,
                 'shelf':shelf, 'peak_props':peak_props, 'peak_props_orig':peak_props_orig, 'peak_ldrs_orig': peak_ldrs_orig}
    if(len(peak_ldrs)==2):
        out={'peak1': peak_ldrs[0], 'peak2': peak_ldrs[1],  'peak1_height':peak_height1, 'peak2_height':peak_height2,
             'shelf':shelf, 'peak_props':peak_props, 'peak_props_orig':peak_props_orig, 'peak_ldrs_orig': peak_ldrs_orig}
    if(len(peak_ldrs)>2):
        print("More than two peaks returned. This shouldn't happen -- must be a bug in the code")
        out={'peak1': peak_ldrs[0], 'peak2': peak_ldrs[1],  'peak1_height':peak_height1, 'peak2_height':peak_height2,
             'shelf':shelf, 'peak_props':peak_props, 'peak_props_orig':peak_props_orig, 'peak_ldrs_orig': peak_ldrs_orig}
    return(out)

def get_ldr_cutoff_mixture(logint, peak_ldrs, show=True, mean_tol=0.4, silent=False):
    #print("get_ldr_cutoff_mixture")
    logint = logint.copy()
    peak_ldrs=peak_ldrs.copy()
    ### mixture model
    X = np.array(logint).reshape(-1,1)
    X = torch.tensor(X).float()
    try:
        #print("model fitting")
        #print(peak_ldrs[0])
        #print(peak_ldrs[1])
        m1 = torch.tensor(peak_ldrs[0])
        m2 = torch.tensor(peak_ldrs[1])
        #m1.frozen=True
        #m2.frozen=True
        d1 = Normal(means=[m1], frozen = False)
        d2 = Normal(means=[m2], frozen = False)
        #d1.frozen=torch.tensor(True)
        #d2.frozen=torch.tensor(True)
        d3 = [d1,d2]
        priors = np.empty((len(X), 2))
        for i in range(len(priors)):
            if X[i][0] < peak_ldrs[0]+0.5:
                priors[i][0] = 1
                priors[i][1] = 0
            elif X[i][0] > peak_ldrs[1]-0.5:
                priors[i][0] = 0
                priors[i][1] = 1
            else:
                priors[i][0] = 0.5
                priors[i][1] = 0.5
        model = GeneralMixtureModel(d3, verbose=False, frozen=False, tol=0.001, max_iter=100, inertia=0.9).fit(X, priors=priors)
    except:
        if not silent: print("mixture model failed")
        return(float('nan'))
    try:
        ldr_cutoff = get_mixture_cutoff(model, silent=True)
    except:
        if not silent: print("get mixture cutoff failed")
        return(float('nan'))
    if show:
        try:
            x = np.arange(np.min(logint), np.max(logint), 0.1)
            y1 = model.distributions[0].probability(x.reshape(-1, 1))
            y2 = model.distributions[1].probability(x.reshape(-1, 1))
            y3 = model.probability(x.reshape(-1, 1))
            #fig, ax = plt.subplots()
            plt.figure(figsize=(6, 3))
            plt.hist(X[:,0], density=True, bins=30)
            plt.plot(x, y1, color = "green", label="Normal1")
            plt.axvline(peak_ldrs[0], color="green", label="live peak")
            plt.plot(x, y2, color = "red", label="Normal2")
            plt.axvline(peak_ldrs[1], color="red", label="dead peak")
            plt.plot(x, y3, color = "purple", label="Mixture")
            plt.axvline(ldr_cutoff, color="orange", label="LDR cutoff")
            plt.legend(loc=(1.05, 0.4))
            plt.tight_layout()
            plt.show()
        except:
            if not silent: print("plotting mixture model failed")
            return(float('nan'))
    mean1 = model.distributions[0].means.item()
    mean2 = model.distributions[1].means.item()
    check1 = abs(mean1-peak_ldrs[0]) < mean_tol
    check2 = abs(mean2-peak_ldrs[1]) < mean_tol
    check3 = ldr_cutoff > peak_ldrs[0] and ldr_cutoff < peak_ldrs[1]
    if(check1 and check2 and check3):
        return(ldr_cutoff)
    else:
        if not silent: print("mixture model fitting failed")
        return(float('nan'))

def get_ldr_cutoff_valley(x,y, peak_ldrs, silent=False):
    x=x.copy()
    y=y.copy()
    peak_ldrs=peak_ldrs.copy()
    ### find valley in between two most prominent peaks
    x_sub = [val for val in x if val < peak_ldrs[1] and val > peak_ldrs[0] ]
    y_sub = [val for val in y if val < peak_ldrs[1] and val > peak_ldrs[0] ]
    y_sub_neg = [-val for val in y_sub]
    valley_locs, valley_props = find_peaks(y_sub_neg, height=float('-Inf'), prominence=0)
    if len(valley_locs) > 0:
        ### get "peak" with maximum height -- all "peaks" will have negative height, so this will give the lowest valley
        valley_ldrs = [ y_sub[loc] for loc in valley_locs ]
        ldr_cutoff = np.max(valley_ldrs)
        return(ldr_cutoff)
    else:
        if not silent: print("valley method failed")
        return(float('nan'))

def get_ldrgates_new(ldrint, smoothing=1, show=True, first_peak_min=0.5,
                     min_prominence=0, min_peak_height=0.02, min_peak_distance=0.5,
                     single_peak_cutoff=3,
                     mixture_backup_method="valley", silent=True,
                    return_peaks_only=False):
    ldrint = ldrint[ldrint > 0]
    logint = np.log10(ldrint)
    logint = logint[ [not x for x in np.isnan(logint)] ]
    logint = logint[ [not x for x in np.isinf(logint)] ]
    logint = logint.copy()
    if show:
        fig, ax = plt.subplots()
        x, y = sns.kdeplot(logint, ax=ax, bw_adjust=smoothing).get_lines()[0].get_data()
        plt.close()
    else:
        kde = sns.kdeplot(logint, bw_adjust=smoothing)
        x, y = kde.get_lines()[0].get_data()
        plt.close()
    alive_dead_peaks = get_peaks_ldr(x.copy(), y.copy(), smoothing=smoothing, first_peak_min=first_peak_min,
                                    min_prominence=min_prominence, min_peak_height=min_peak_height,
                                    min_peak_distance=min_peak_distance, single_peak_cutoff=single_peak_cutoff,
                                    silent=silent)
    if return_peaks_only: return(alive_dead_peaks)
    peak_ldrs = [alive_dead_peaks['peak1'], alive_dead_peaks['peak2']]
    peak_ldrs = [x for x in peak_ldrs if not np.isnan(x)]
    method_used="none"
    if show:
        fig, ax = plt.subplots()
        x, y = sns.kdeplot(logint, ax=ax, bw_adjust=smoothing).get_lines()[0].get_data()
        #plt.show()
    else:
        kde = sns.kdeplot(logint, bw_adjust=smoothing)
        x, y = kde.get_lines()[0].get_data()
        plt.close()
    if len(peak_ldrs) == 1:
        if not silent: print("one peak. almost all cells alive?")
        ### assume almost all cells are alive
        if peak_ldrs[0] < single_peak_cutoff:
            if not silent: print("one peak. almost all cells alive?")
            ldr_cutoff = np.quantile(logint, 0.99)
            ldr_cutoff_mixture=ldr_cutoff
            ldr_cutoff_valley=ldr_cutoff
            ldr_cutoff_middle=ldr_cutoff
        else:
            if not silent: print("one peak. almost all cells dead?")
            ldr_cutoff = np.quantile(logint, 0.01)
            ldr_cutoff_mixture=ldr_cutoff
            ldr_cutoff_valley=ldr_cutoff
            ldr_cutoff_middle=ldr_cutoff
    elif len(peak_ldrs) > 1:
        if return_peaks_only: return({'peak_ldrs':peak_ldrs, 'peak_props':peak_props})
        ###### Note: write up this section as a new function: "get_ldr_cutoff"
        ldr_cutoff_mixture = get_ldr_cutoff_mixture(logint.copy(), peak_ldrs, show=show, silent=silent)
        ldr_cutoff_valley = get_ldr_cutoff_valley(x,y, peak_ldrs, silent=silent)
        ldr_cutoff_middle = (peak_ldrs[0] + peak_ldrs[1])/2
        if not np.isnan(ldr_cutoff_mixture):
            method_used = "mixture"
            ldr_cutoff = ldr_cutoff_mixture
        elif not np.isnan(ldr_cutoff_valley):
            method_used = "valley"
            ldr_cutoff = ldr_cutoff_valley
        else:
            method_used = "middle"
            ldr_cutoff = ldr_cutoff_middle
    out = alive_dead_peaks
    out['method_used'] = method_used
    out['ldr_cutoff'] = ldr_cutoff
    out['ldr_gates'] = np.array([-np.inf, ldr_cutoff])
    out['ldr_cutoff_mixture'] = ldr_cutoff_mixture
    out['ldr_cutoff_valley'] = ldr_cutoff_valley
    out['ldr_cutoff_middle'] = ldr_cutoff_middle
    #ldr_lims = np.array([x.min(), x.max()])
    #return(ldr_gates, ldr_lims)
    if not silent: print(out)
    return(out)

def get_ldrgates_new_well(barcode, well, smoothing=1.1, show=True, first_peak_min=0.5,
                          min_prominence=0, min_peak_height=0.02, min_peak_distance=0.5,
                          single_peak_cutoff=3,
                          silent=True, return_peaks_only=False):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    df = read_and_rename_well_data(barcode, well, silent=True)
    ldrint = df['ldr'].copy()
    results = get_ldrgates_new(ldrint, smoothing=smoothing, show=show, silent=silent,
                               min_prominence=min_prominence, min_peak_height=min_peak_height,
                               min_peak_distance=min_peak_distance, single_peak_cutoff=single_peak_cutoff,
                               return_peaks_only=return_peaks_only)
    return(results)


def get_peaks_cell_line(barcode, cell_line, method="mixture", smoothing=1, df_only=True, return_peaks_only=False):
    wells = get_wells(barcode, cell_line)
    peak_list=[]
    df_list = []
    for well in wells:
        print(well)
        pdict = get_ldrgates_new_well(barcode, well, show=True, return_peaks_only=return_peaks_only,
                                      smoothing=smoothing, silent=True)
        pdict['barcode']=barcode
        pdict['cell_line']=cell_line
        pdict['well']=well
        if return_peaks_only:
            cols = ['barcode', 'cell_line', 'well', 'peak1', 'peak2', 
                    'peak1_height', 'peak2_height','shelf']
        else:
            cols = ['barcode', 'cell_line', 'well', 'ldr_cutoff', 'peak1', 'peak2', 
                    'peak1_height', 'peak2_height','shelf','method_used', 'ldr_cutoff_mixture', 'ldr_cutoff_valley', 'ldr_cutoff_middle']
        dd = {k:v for (k,v) in pdict.items() if k in cols}
        df = pd.DataFrame(data=dd, index=[0])
        df = df[cols]
        df_list.append(df)
        peak_list.append(pdict['peak_props'])
    df_full = pd.concat(df_list)
    df_full = df_full.reset_index()
    if df_only:
        return(df_full)
    else:
        ### return the peak lists later
        return({'df':df_full, 'peak_props':peak_list})

def get_peaks_barcode(barcodes, smoothing=1, df_only=True, return_peaks_only=False):
    peak_list=[]
    df_list = []
    for barcode in barcodes:
        cell_lines = get_cell_lines_on_plate(barcode)
        for cell_line in cell_lines:
            #print(cell_line)
            wells = get_wells(barcode, cell_line)
            for well in wells:
                pdict = get_ldrgates_new_well(barcode, well, show=False, return_peaks_only=return_peaks_only,
                                      smoothing=smoothing, silent=True)
                pdict['barcode']=barcode
                pdict['cell_line']=cell_line
                pdict['well']=well
                if return_peaks_only:
                    cols = ['barcode', 'cell_line', 'well', 'peak1', 'peak2', 
                            'peak1_height', 'peak2_height','shelf']
                else:
                    cols = ['barcode', 'cell_line', 'well', 'ldr_cutoff', 'peak1', 'peak2', 
                            'peak1_height', 'peak2_height','shelf','method_used', 'ldr_cutoff_mixture', 'ldr_cutoff_valley', 'ldr_cutoff_middle']
                dd = {k:v for (k,v) in pdict.items() if k in cols}
                df = pd.DataFrame(data=dd, index=[0])
                df = df[cols]
                df_list.append(df)
                peak_list.append(pdict['peak_props'])
    df_full = pd.concat(df_list)
    df_full = df_full.reset_index()
    if df_only:
        return(df_full)
    else:
        ### return the peak lists later
        return({'df':df_full, 'peak_props':peak_list})

def get_all_peaks(smoothing=1):
    load_well_metadata('meta')
    barcodes = meta['barcode'].unique()
    out_file = "peak_locations_new.csv"
    if os.path.exists(out_file):
        df = pd.read_csv(out_file)
        done_barcodes = list(df['barcode'])
        barcodes = [bc for bc in barcodes if not (bc in done_barcodes)]
    barcodes = sorted(barcodes)
    for i in range(len(barcodes)):
        barcode = barcodes[i]
        print(barcode)
        res = get_peaks_barcode([barcode], smoothing=smoothing)
        if i==0 and not ('done_barcodes' in locals() or 'done_barcodes' in globals()):
            res.drop(columns=['index']).to_csv(out_file, header=True, index=False)
        else:
            res.drop(columns=['index']).to_csv(out_file, mode='a', header=False, index=False)
    return(res)

## steps
# 1) identify ends of left and right shelves
# 2) if main peak < 3 (default) assume it's the live peak, otherwise assume it's the dead peak
# 3) identify all x intervals that are:
#        1) inbetween those two shelves
#        2) at least min_distance=1 from first peak
#        3) to the right of the peak (if main peak<3), or to the left of the peak  (if main peak < 3)
# 4) identify the point closest to zero
# 3) report the point as the second peak as long as it's close enough to 0

### input x and y from the sns.kdeplot function
def find_side_shelves(x, y, tol=0.03):
    x = x.copy()
    y = y.copy()
    dx = np.diff(x)[0]
    grad = np.gradient(y, dx)  ## first derivative
    #grad2 = np.gradient(grad, dx)  ## second derivative
    for i in range(len(grad)):
        if abs(grad[i]) > tol:
            end_left_shelf = x[i]
            break
    for i in reversed(range(len(grad))):
        if abs(grad[i]) > tol:
            end_right_shelf = x[i]
            break
    return({'end_left':end_left_shelf, 'end_right':end_right_shelf})

def find_shelf_old(x, y, main_peak_ldr, min_peak_distance=1, single_peak_cutoff=3, tol=0.03):
    sides = find_side_shelves(x,y,tol)
    dx = np.diff(x)[0]
    grad = np.gradient(y, dx)  ## first derivative
    #grad2 = np.gradient(grad, dx)  ## second derivative
    ### remove side shelves
    ind = [list(x).index(val) for val in x if val > sides['end_left'] and val < sides['end_right']]
    x = [x[i] for i in ind]
    y = [y[i] for i in ind]
    grad = [grad[i] for i in ind]
    ### remove values close to the main peak (and to the wrong side of peak)
    if main_peak_ldr < single_peak_cutoff: ## assume it's alive peak
        ind = [list(x).index(val) for val in x if val > main_peak_ldr + min_peak_distance]
        x = [x[i] for i in ind]
        y = [y[i] for i in ind]
        grad = [grad[i] for i in ind]
    else: ### assume it's the dead peak
        ind = [list(x).index(val) for val in x if val < main_peak_ldr - min_peak_distance]
        x = [x[i] for i in ind]
        y = [y[i] for i in ind]
        grad = [grad[i] for i in ind]
    ind_shelf = np.argmin(np.abs(grad))
    #print(grad[ind_shelf])
    if np.abs(grad[ind_shelf]) < tol:
        #print('found')
        return(x[ind_shelf])
    else:
        #print("not found")
        return(None)

def find_shelf(x, y, main_peak_ldr, min_peak_distance=1, single_peak_cutoff=3, tol=0.03, 
               min_peak_height=0.02, first_peak_min=0.5, max_slope=0.3):
    x = x.copy()
    y = y.copy()
    sides = find_side_shelves(x,y,tol)
    dx = np.diff(x)[0]
    grad = np.gradient(y, dx)  ## first derivative
    #return(grad)
    #grad2 = np.gradient(grad, dx)  ## second derivative
    ### remove values close to the main peak (and to the wrong side of peak)
    if main_peak_ldr < single_peak_cutoff: ## assume it's alive peak
        #print('alive')
        #print(main_peak_ldr + min_peak_distance)
        #print(sides['end_right'])
        ind = [list(x).index(val) for val in x if val > main_peak_ldr + min_peak_distance and val < sides['end_right']]
        x = [x[i] for i in ind]
        y = [y[i] for i in ind]
        grad = [grad[i] for i in ind]
        ### find peaks close to zero
        peak_locs, peak_props = find_peaks(grad, height=float('-Inf'), prominence=0)
    else: ### assume it's the dead peak
        #print('dead')
        #print(main_peak_ldr - min_peak_distance)
        #print(sides['end_left'])
        ind = [list(x).index(val) for val in x if val < main_peak_ldr - min_peak_distance and val > sides['end_left']]
        x = [x[i] for i in ind]
        y = [y[i] for i in ind]
        grad = [grad[i] for i in ind]
        neg_grad = [-x for x in grad]
        ### find valleys close to zero
        peak_locs, peak_props = find_peaks(neg_grad, height=float('-Inf'), prominence=0)
    #return([x[i] for i in peak_locs])
    if len(peak_locs)>0:
        max_ind = np.argmin(abs(np.array(peak_props['peak_heights'])))
        ### choose the peak with height closed to zero
        shelf_loc = peak_locs[max_ind]
        ldr_shelf = x[shelf_loc] ### return x coordinate (ldr intensity)
        shelf_height = y[shelf_loc]
        shelf_slope = grad[shelf_loc]
        if ldr_shelf > first_peak_min and shelf_height>min_peak_height and abs(shelf_slope) < max_slope:
            out = {'ldr':ldr_shelf, 'height':shelf_height, 'slope': shelf_slope}
            #print(out)
            return(out)
        else:
            return(None)
    else:
        return(None)

def plot_density_derivative(barcode, well):
    logint = get_logldrint(barcode, well)
    fig, ax = plt.subplots()
    x, y = sns.kdeplot(logint, ax=ax, bw_adjust=1).get_lines()[0].get_data()
    plt.close()
    #sides = find_side_shelves(x,y,tol)
    dx = np.diff(x)[0]
    grad = np.gradient(y, dx)  ## first derivative
    gg1 = ggplot()+geom_line(aes(x=x, y=grad))
    print(gg1)
    out = {'x': x,'grad': grad}
    neg_grad = [-x for x in grad]
    peak_locs, peak_props = find_peaks(-grad, height=float('-Inf'), prominence=0)
    out = {'peak_locs': peak_locs, 'peak_props': peak_props}
    return(out)

def get_skewness(barcode, well):
    logint = get_logldrint(barcode, well, remove_na=True)
    fig, ax = plt.subplots()
    x, y = sns.kdeplot(logint, ax=ax, bw_adjust=1).get_lines()[0].get_data()
    #plt.close()
    skewness = skew(logint)
    test_res = skewtest(logint)
    return({'skewness':skewness, 'statistic': test_res[0], 'pvalue': test_res[1]})

### function to get all of the ldrint values and gate them as one distribution
def get_density_together_single_dna(barcode, cell_line, plot_wells=True):
    if not isinstance(barcode, str):
        print(barcode)
        print(cell_line)
        print("barcode needs to be a string")
        return(None)
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    wells = get_wells(barcode, cell_line)
    df_list=[]
    for well in wells:
        df = read_and_rename_well_data(barcode, well, silent = True)
        df['well'] = well
        df_list.append(df)
    df = pd.concat(df_list)
    #res = get_ldrgates_new(df['dna'].copy(), show=False, silent=True)
    fig, ax = plt.subplots()
    if plot_wells:
        for well in wells:
            #print(well)
            df_sub = df.query("well == '" + well + "'")
            dna = df_sub['dna'].copy()
            #ldrint = ldrint[ldrint > 0]
            logdna = np.log10(dna)
            #logint = logint[ [not x for x in np.isnan(logint)] ]
            #logint = logint[ [not x for x in np.isinf(logint)] ]
            sns.kdeplot(logdna, ax=ax, alpha=0.25).set_title(cell_line)
    dna_all = df['dna'].copy()
    #ldrint = ldrint[ldrint > 0]
    #logint = np.log10(ldrint)
    #logint = logint[ [not x for x in np.isnan(logint)] ]
    #logint = logint[ [not x for x in np.isinf(logint)] ]
    sns.kdeplot(np.log10(dna_all), ax=ax, alpha=1, color="black", linewidth=2).set_title(cell_line)
    #plt.axvline(x=res['peak1'], ls = "--", color = "green")
    #plt.axvline(x=res['ldr_cutoff'], ls = "--", color = "orange")
    #plt.axvline(x=res['peak2'], ls = "--", color = "red")
    plt.show()
    return(res)

    
### function to get all of the ldrint values and gate them as one distribution
def get_density_together_single(barcode, cell_line, plot_wells=True):
    if not isinstance(barcode, str):
        print(barcode)
        print(cell_line)
        print("barcode needs to be a string")
        return(None)
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    wells = get_wells(barcode, cell_line)
    df_list=[]
    for well in wells:
        df = read_well_data(barcode, well)
        df['well'] = well
        df_list.append(df)
    df = pd.concat(df_list)
    res = get_ldrgates_new(df['ldrint'].copy(), show=False, silent=True)
    fig, ax = plt.subplots()
    if plot_wells:
        for well in wells:
            df_sub = df.query("well == '" + well + "'")
            ldrint = df_sub['ldrint'].copy()
            ldrint = ldrint[ldrint > 0]
            logint = np.log10(ldrint)
            logint = logint[ [not x for x in np.isnan(logint)] ]
            logint = logint[ [not x for x in np.isinf(logint)] ]
            sns.kdeplot(logint, ax=ax, alpha=0.25).set_title(cell_line)
    ldrint = df['ldrint'].copy()
    ldrint = ldrint[ldrint > 0]
    logint = np.log10(ldrint)
    logint = logint[ [not x for x in np.isnan(logint)] ]
    logint = logint[ [not x for x in np.isinf(logint)] ]
    sns.kdeplot(logint, ax=ax, alpha=1, color="black", linewidth=2).set_title(cell_line)
    plt.axvline(x=res['peak1'], ls = "--", color = "green")
    plt.axvline(x=res['ldr_cutoff'], ls = "--", color = "orange")
    plt.axvline(x=res['peak2'], ls = "--", color = "red")
    plt.show()
    return(res)

def get_density_together_barcode(barcode, plot_wells=True, filename = "pdf_test.pdf", no_pdf=False):
    if not isinstance(barcode, str):
        print(barcode)
        print(cell_line)
        print("barcode needs to be a string")
        return(None)
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    #fig, ax = plt.subplots()
    cell_lines = get_cell_lines_on_plate(barcode)
    res_list = []
    pdf_full = filename
    nb_plots = len(cell_lines)
    if nb_plots != 6: print(barcode + ": " + str(nb_plots) + " cell lines")
    if nb_plots == 0: return(None)
    ncols = 2
    nrows = int(np.ceil(nb_plots/2))
    res_list = []
    df_full_list = []
    for cell_line in cell_lines:
        #print(cell_line)
        wells = get_wells(barcode, cell_line)
        df_list=[]
        for well in wells:
            df = read_well_data(barcode, well)
            df['well'] = well
            df_list.append(df)
        df = pd.concat(df_list)
        df_full_list.append(df)
        pdict = get_ldrgates_new(df['ldrint'].copy(), show=False, silent=True)
        pdict['barcode']=barcode
        pdict['cell_line']=cell_line
        cols = ['barcode', 'cell_line', 'ldr_cutoff', 'peak1', 'peak2', 
                    'peak1_height', 'peak2_height','shelf','method_used', 'ldr_cutoff_mixture', 'ldr_cutoff_valley', 'ldr_cutoff_middle']
        dd = {k:v for (k,v) in pdict.items() if k in cols}
        res_df = pd.DataFrame(data=dd, index=[0])
        res_df = res_df[cols]
        res_list.append(res_df)
    res_df_full = pd.concat(res_list)
    res_df_full = res_df_full.reset_index(drop=True)
    if no_pdf:
        return(res_df_full)
    pdf_pages = PdfPages(pdf_full)
    figs, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(9, 4*nrows),
                        layout="constrained", sharex= "all")
    figs.suptitle(barcode, fontsize=16)
    for i, ax in zip(range(len(cell_lines)), axes.ravel()):
        cell_line = cell_lines[i]
        df = df_full_list[i]
        res = res_list[i]
        wells = get_wells(barcode, cell_line)
        #fig, ax = plt.subplots()
        #fig = plt.figure()
        if plot_wells:
            for well in wells:
                df_sub = df.query("well == '" + well + "'")
                ldrint = df_sub['ldrint'].copy()
                ldrint = ldrint[ldrint > 0]
                logint = np.log10(ldrint)
                logint = logint[ [not x for x in np.isnan(logint)] ]
                logint = logint[ [not x for x in np.isinf(logint)] ]
                sns.kdeplot(logint, ax=ax, alpha=0.25).set_title(cell_line)
        ldrint = df['ldrint'].copy()
        ldrint = ldrint[ldrint > 0]
        logint = np.log10(ldrint)
        logint = logint[ [not x for x in np.isnan(logint)] ]
        logint = logint[ [not x for x in np.isinf(logint)] ]
        sns.kdeplot(logint, ax=ax, alpha=1, color="black", linewidth=2).set_title(cell_line)
        ax.axvline(x=res['peak1'], ls = "--", color = "green")
        ax.axvline(x=res['ldr_cutoff'], ls = "--", color = "orange")
        ax.axvline(x=res['peak2'], ls = "--", color = "red")
        #figures.append(fig)
        #plt.show()
    #plt.show()
    #for i, ax in enumerate(axes.flat):
    #    if i < len(figures):
    #        ax.imshow(figures[i])  # Adjust as needed
    #        #ax.set_title("Figure {}".format(i + 1))  # Adjust title
    plt.tight_layout()
    pdf_pages.savefig(figs)
    plt.close('all')
    pdf_pages.close()
    return(res_df_full)

def plot_density_together_all():
    load_well_metadata('meta')
    barcodes = meta['barcode'].unique()
    barcodes.sort()
    for barcode in barcodes:
        ff = os.path.join("kde_plots_joint", barcode + ".pdf")
        print(barcode)
        if not os.path.isfile(ff):
            get_density_together_barcode(barcode, filename = ff)

def get_ldr_cutoffs_together_all(out_file="new_gating_together.csv"):
    load_well_metadata('meta')
    barcodes = meta['barcode'].unique()
    if os.path.exists(out_file):
        df_old = pd.read_csv(out_file)
        done_barcodes = list(df_old['barcode'])
        barcodes = [bc for bc in barcodes if not (bc in done_barcodes)]
    barcodes = sorted(barcodes)
    df_list=[]
    for i in range(len(barcodes)):
        barcode = barcodes[i]
        print(barcode)
        df = get_density_together_barcode(barcode, no_pdf = True)
        df_list.append(df)
        if i==0 and not ('done_barcodes' in locals() or 'done_barcodes' in globals()):
            df.to_csv(out_file, header=True, index=False)
        else:
            df.to_csv(out_file, mode='a', header=False, index=False)
            
    df = pd.concat(df_list)
    df = df.reset_index(drop=True)
    return(df)

def barcode_from_number(number):
    if not "bc_to_num" in globals(): define_bc_to_num('bc_to_num')
    return(bc_to_num[str(number)])

def define_bc_to_num(name='bc_to_num'):
    bc_to_num = {
        '12': '201117_combo_12','13': '201117_combo_13','14': '201117_combo_14','15': '201117_combo_15',
        '16': '201117_combo_16','17': '201117_combo_17','18': '201117_combo_18','19': '201117_combo_19',
        '20': '201117_combo_20','21': '201117_combo_21','25': '201117_combo_25','26': '201117_combo_26',
        '27': '201117_combo_27','28': '201117_combo_28','29': '201117_combo_29','30': '201117_combo_30',
        '31': '201117_combo_31','32': '201117_combo_32','33': '201117_combo_33','34': '201117_combo_34',
        '38': '210219_combo_38','39': '210219_combo_39','40': '210219_combo_40','41': '210219_combo_41',
        '42': '210219_combo_42','43': '210219_combo_43','44': '210219_combo_44','45': '210219_combo_45',
        '46': '210219_combo_46','47': '210219_combo_47','48': '210226_combo_48','49': '210226_combo_49',
        '50': '210226_combo_50','51': '210226_combo_51','52': '210226_combo_52','53': '210226_combo_53',
        '54': '210226_combo_54','55': '210226_combo_55','56': '210226_combo_56','57': '210226_combo_57',
        '58': '210302_combo_58','59': '210302_combo_59','60': '210302_combo_60','61': '210302_combo_61',
        '62': '210406_combo_62','69': '210406_combo_69','70': '210406_combo_70','71': '210406_combo_71',
        '72': '210406_combo_72','73': '210406_combo_73','74': '210406_combo_74','75': '210406_combo_75',
        '76': '210406_combo_76','77': '210406_combo_77','78': '210423_combo_78','79': '210423_combo_79',
        '80': '210423_combo_80','81': '210423_combo_81','82': '210423_combo_82','83': '210423_combo_83',
        '84': '210423_combo_84','85': '210423_combo_85','86': '210423_combo_86','87': '210423_combo_87',
        '88': '210518_combo_88','89': '210518_combo_89','90': '210518_combo_90','91': '210518_combo_91',
        '92': '210518_combo_92','93': '210518_combo_93','94': '210518_combo_94','95': '210518_combo_95',
        '96': '210518_combo_96','97': '210518_combo_97','100': '210521_combo_100','101': '210521_combo_101',
        '102': '210521_combo_102','103': '210521_combo_103','104': '210521_combo_104','105': '210521_combo_105',
        '106': '210521_combo_106','107': '210521_combo_107','98': '210521_combo_98','99': '210521_combo_99',
        '108': '210611_combo_108','109': '210611_combo_109','110': '210611_combo_110','111': '210611_combo_111',
        '112': '210611_combo_112','113': '210611_combo_113','114': '210611_combo_114','115': '210611_combo_115',
        '116': '210611_combo_116','117': '210611_combo_117','118': '210727_combo_118','119': '210727_combo_119',
        '120': '210727_combo_120','121': '210727_combo_121','122': '210727_combo_122','123': '210727_combo_123',
        '124': '210727_combo_124','125': '210727_combo_125','126': '210727_combo_126','127': '210727_combo_127',
        '128': '210727_combo_128','129': '210727_combo_129','130': '210727_combo_130','131': '210730_combo_131',
        '132': '210730_combo_132','133': '210730_combo_133','134': '210730_combo_134','135': '210730_combo_135',
        '136': '210730_combo_136','137': '210730_combo_137','138': '210730_combo_138','139': '210730_combo_139',
        '140': '210730_combo_140','141': '210730_combo_141','142': '210730_combo_142','143': '210730_combo_143',
        '144': '210806_combo_144','145': '210806_combo_145','146': '210806_combo_146','147': '210806_combo_147',
        '148': '210806_combo_148','149': '210806_combo_149','150': '210806_combo_150','151': '210806_combo_151',
        '152': '210806_combo_152','155': '210806_combo_155','156': '210806_combo_156','157': '211005_combo_157',
        '158': '211005_combo_158','159': '211005_combo_159','160': '211005_combo_160','161': '211005_combo_161',
        '162': '211005_combo_162','163': '211005_combo_163','164': '211005_combo_164','165': '211005_combo_165',
        '166': '211005_combo_166','167': '211015_combo_167','168': '211015_combo_168','169': '211015_combo_169',
        '170': '211015_combo_170','171': '211015_combo_171','172': '211015_combo_172','173': '211015_combo_173',
        '174': '211015_combo_174','175': '211015_combo_175','176': '211015_combo_176','177': '211029_combo_177',
        '178': '211029_combo_178','179': '211029_combo_179','180': '211029_combo_180','181': '211029_combo_181',
        '182': '211029_combo_182','183': '211029_combo_183','184': '211029_combo_184','185': '211029_combo_185',
        '186': '211029_combo_186'
    }
    globals()[name] = bc_to_num
    return(None)


def get_ldrgates_joint(ldrint_list, wells, select = "median", window_size = 0.2):
    """Joint gating of many wells, based on ldr intensities
    This function gates each well separately and then 

    Parameters
    ----------
    ldrint_list : list of lists
        Each item in the list should be an ldrint feature (a list itself) across all cells in a well
    wells: list
        a list of the names of each well: e.g. "M07", "C03", etc. -- each name must be unique
    select: a string
        default is "median" -- construct a window around the median ldr_cutoff and snap all cutoffs to this window.
    window_size: float
        default is 0.2 -- construct a window of +/- 0.2 around the median ldr_cutoff

    Returns
    -------
    results : dictionary
        The dictionary contains a dataframe of results and a list of peak properties for each well.
        The keys are:
            df: pandas data frame
                Each row contains the live/dead peak locations and ldr_cutoff location, etc. for one well.
            peak_props: a list
                Each item is a dictionary containing peak properties for each well.
    """
    ### write a check that no well names are duplicated: return error if any are
    
    df_list = []
    props_dict = {}
    for i in range(len(ldrint_list)):
        ldrint = ldrint_list[i]
        well = wells[i]
        pdict = get_ldrgates_new(ldrint.copy(), show=False, silent=True)
        pdict['well'] = well
        #pdict['barcode']=barcode
        #pdict['cell_line']=cell_line
        #cols = ['barcode', 'cell_line', 'ldr_cutoff', 'peak1', 'peak2', 
        #            'peak1_height', 'peak2_height','shelf','method_used', 'ldr_cutoff_mixture', 
        #        'ldr_cutoff_valley', 'ldr_cutoff_middle']
        cols = ['well','ldr_cutoff', 'peak1', 'peak2', 
                    'peak1_height', 'peak2_height','shelf','method_used', 'ldr_cutoff_mixture', 
                'ldr_cutoff_valley', 'ldr_cutoff_middle']
        dd = {k:v for (k,v) in pdict.items() if k in cols}
        pdict2 = {k:v for (k,v) in pdict.items() if not k in cols}
        res_df = pd.DataFrame(data=dd, index=[0])
        res_df = res_df[cols]
        df_list.append(res_df)
        props_dict[well] = pdict2
    res_df_full = pd.concat(df_list)
    res_df_full = res_df_full.reset_index(drop=True)
    ldr_cutoff_median = np.nanmedian(res_df_full.ldr_cutoff)
    ldr_cutoff_min = ldr_cutoff_median - window_size
    ldr_cutoff_max = ldr_cutoff_median + window_size
    res_df_full['ldr_cutoff_median'] = ldr_cutoff_median
    res_df_full['ldr_cutoff_min'] = ldr_cutoff_min
    res_df_full['ldr_cutoff_max'] = ldr_cutoff_max
    res_df_full['ldr_cutoff_final'] = np.clip(res_df_full['ldr_cutoff'], ldr_cutoff_min, ldr_cutoff_max)
    results = {'df': res_df_full, 'peak_props': props_dict}
    return(results)

def get_counts_new(barcode, cell_lines=None, plot_wells=True, select = "median", window_size = 0.2, 
                   n_wells = None):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    if cell_lines is None:
        cell_lines = get_cell_lines_on_plate(barcode)
    else:
        if isinstance(cell_lines, str):
            cell_lines = [cell_lines]
    res_dict_full = {}
    res_df_list = []
    ### get ldr cutoffs for each well in each cell line on the plate
    df_ldr_full_list = []
    for cl in cell_lines:
        #print(cl)
        wells = get_wells(barcode, cl)
        if n_wells is not None:
            wells = wells[0:n_wells] ### try first n wells for testing
        ldrint_list = []
        df_list=[]
        #print(wells)
        dna_colnames = []
        for well in wells:
            #print(well)
            #df = read_well_data(barcode, well)
            df = read_and_rename_well_data(barcode, well, silent = True)
            df['well'] = well
            dna_colnames.append(df['dna_colname'][0])
            ldrint_list.append( df['ldr'].copy() )
            df_list.append(df)
        df = pd.concat(df_list)
        df_ldr_full_list.append(df)
        ### get ldr cutoffs for all wells for one cell line
        ### (gating wells separately, then restrict to a window around the median ldr cutoff)
        results = get_ldrgates_joint(ldrint_list = ldrint_list, wells = wells, 
                                     select = select, window_size = window_size)
        res_df = results['df']
        res_df['barcode'] = barcode
        res_df['cell_line'] = cl
        res_df['dna_colname'] = dna_colnames
        res_df_list.append(res_df)
        res_dict = results['peak_props']
        res_dict_full |= res_dict
        #out_dict = {(barcode + k):v for (k,v) in res_dict.items()}
        #joint_res = get_ldrgates_new(df['ldrint'].copy(), show=False, silent=True)
        if plot_wells: fig, ax = plt.subplots()
        
        for well in wells:
            ### get ldrint for each individual well
            df_sub = df.query("well == '" + well + "'")
            ldrint = df_sub['ldr'].copy()
            ldrint = ldrint[ldrint > 0]
            logint = np.log10(ldrint)
            logint = logint[ [not x for x in np.isnan(logint)] ]
            logint = logint[ [not x for x in np.isinf(logint)] ]
            if plot_wells:
                sns.kdeplot(logint, ax=ax, alpha=0.25).set_title(cl)
        ### get ldrint for all wells together
        ldrint = df['ldr'].copy()
        ldrint = ldrint[ldrint > 0]
        logint = np.log10(ldrint)
        logint = logint[ [not x for x in np.isnan(logint)] ]
        logint = logint[ [not x for x in np.isinf(logint)] ]
        if plot_wells:
            sns.kdeplot(logint, ax=ax, alpha=1, color="black", linewidth=2).set_title(cl)
            plt.axvline(x=np.nanmedian(res_df['peak1']), ls = "--", color = "green")
            plt.axvline(x=res_df['ldr_cutoff_median'][0], ls = "-", color = "orange")
            plt.axvline(x=res_df['ldr_cutoff_min'][0], ls = "--", color = "orange")
            plt.axvline(x=res_df['ldr_cutoff_max'][0], ls = "--", color = "orange")
            plt.axvline(x=np.nanmedian(res_df['peak2']), ls = "--", color = "red")
            plt.show()
    out_df = pd.concat(res_df_list)
    out_df = out_df.reset_index(drop=True)
    out_df.rename(columns = {'ldr_cutoff': 'ldr_cutoff_orig'}, inplace = True)
    out_df['window_size'] = window_size
    out_df['select'] = select
    out_df['snapped'] = np.where(out_df['ldr_cutoff_orig'] != out_df['ldr_cutoff_final'], True, False)
    cols = ['barcode', 'cell_line', 'well', 'ldr_cutoff_final', 'ldr_cutoff_orig', 'ldr_cutoff_median',
        'window_size', 'select', 'peak1', 'peak2',
       'peak1_height', 'peak2_height', 'shelf', 'method_used', 'ldr_cutoff_mixture',
       'ldr_cutoff_valley', 'ldr_cutoff_middle', 'ldr_cutoff_min', 'ldr_cutoff_max', 'snapped', 'dna_colname']
    out_df = out_df[cols]
    
    #### dna-gating and live-dead calling
    df_ldr_full = pd.concat(df_ldr_full_list)
    df_ldr_full = df_ldr_full.reset_index(drop=True)
    df_counts = []
    for i in range(out_df.shape[0]):
        ldr_cutoff = out_df.ldr_cutoff_final[i]
        well = out_df.well[i]
        df_ldr = df_ldr_full.query("well == '" + well + "'")
        ldrint = df_ldr['ldr'].copy()
        #############
        ldr_gates = np.array([-np.inf, ldr_cutoff])
        ldrint[ldrint < 0] = float('nan')
        logint = np.log10(ldrint)
        logint[np.isnan(logint)] = -10 ### dummy value
        ldr_inner = ((ldr_gates[1] >= logint) & (logint >= ldr_gates[0]))
        # if np.sum(ldr_inner) < 50:
        #     dna = None
        #     dna_gates = None
        # else:
        #     dna = df_ldr['dna'].copy()
        #     dna_gates = dcf_int.get_dna_gating(dna.copy(), ldrint, ldr_gates)
        #     if dna_gates is None:
        #         dna=None
        # dna_colname = list(set(df_ldr['dna_colname']))
        # if len(dna_colname) > 1:
        #     ex = "More than one column name for dna:" + " ".join(dna_colname)
        #     raise Exception(ex)
        # else:
        #     if "Hoechst" in dna_colname: ### if Hoechst, don't do dna gating
        #         dna_gates = None
        dna_gates = None
        dna = None
        cell_fate_dict, outcome = dcf_int.live_dead(ldrint, ldr_gates=ldr_gates, dna=dna, dna_gates= dna_gates)
        live_cols = [s for s in list(cell_fate_dict.keys()) if 'alive' in s]
        dead_cols = [s for s in list(cell_fate_dict.keys()) if 'dead' in s]
        a = 0
        d = 0
        for col in live_cols:
            a += cell_fate_dict[col]
        for col in dead_cols:
            d += cell_fate_dict[col]
        dfs = pd.DataFrame(cell_fate_dict, index=[0])
        dfs.insert(0, 'barcode',  barcode)
        dfs.insert(1, 'well',  well)
        #print('dna_gates')
        #print(dna_gates)
        dfs['cell_count'] = a
        dfs['cell_count__dead'] = d
        dfs['cell_count__total'] = len(ldrint)
        # if dna_gates is not None:
        #     dfs['dna_gate1'] = dna_gates[0]
        #     dfs['dna_gate2'] = dna_gates[1]
        #     dfs['dna_gate3'] = dna_gates[2]
        #     dfs['dna_gate4'] = dna_gates[3]
        # else:
        #     dfs['dna_gate1'] = np.nan
        #     dfs['dna_gate2'] = np.nan
        #     dfs['dna_gate3'] = np.nan
        #     dfs['dna_gate4'] = np.nan
        #dfs['ldr_cutoff'] = ldr_gates[1]
        df_counts.append(dfs)
        #################
    df_counts = pd.concat(df_counts)
    df_counts.reset_index(inplace=True, drop = True)
    out_df = pd.merge(out_df, df_counts, how = "left", on = ["barcode", "well"])

    #### return data frame with gates
    out = {'df': out_df, 'peak_props': res_dict_full} ### note: re-format this dict into a dataframe
    return(out)


def get_counts_all(out_file="lincs_combos_new_counts_2023_08_23", out_dir = "results", 
                   n_wells=None, n_plates=None):
    out_csv = os.path.join(out_dir, out_file + ".csv")
    out_parquet = os.path.join(out_dir, out_file + ".parquet")
    load_well_metadata('meta')
    barcodes = meta['barcode'].unique()
    barcodes = sorted(barcodes)
    if os.path.exists(out_csv):
        df_old = pd.read_csv(out_csv)
        done_barcodes = list(df_old['barcode'])
        barcodes = [bc for bc in barcodes if not (bc in done_barcodes)]
        df_list = [df_old]
    else:
        df_list=[]
    if n_plates is not None:
        rr = range(n_plates)
    else:
        rr = range(len(barcodes))
    for i in rr:
        barcode = barcodes[i]
        print(barcode)
        res = get_counts_new(barcode, plot_wells = False, n_wells = n_wells)
        df = res['df']
        df_list.append(df)
        if i==0 and not ('done_barcodes' in locals().keys() or 'done_barcodes' in globals().keys()):
            df.to_csv(out_csv, header=True, index=False)
        else:
            df.to_csv(out_csv, mode='a', header=False, index=False)
            
    df = pd.concat(df_list)
    df.reset_index(drop=True, inplace=True)
    df.to_parquet(out_parquet)
    return(df)

def kde_dna_plot_wells(barcode, wells, title = "", add_legend=False, smoothing=1, column = "ldr"):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    fig, ax = plt.subplots()
    #ax.set_title(title, fontsize=12)
    for well in wells:
        df = read_and_rename_well_data(barcode, well, silent = True)
        #ldrint = df[column].copy()
        #ldrint[ldrint < 0] = float('nan')
        #logdna = np.log10(ldrint)
        dna = df['dna'].copy()
        logdna = np.log10(dna)
        dna_colname = list(set(df['dna_colname']))[0]
        #x, y = sns.kdeplot(logdna, ax=ax, color = "grey", alpha=0.5).get_lines()[0].get_data()
        if add_legend:
            meta_sub=get_well_meta(barcode=barcode, well=well)
            trt = list(meta_sub.agent1)[0] + " " + list(meta_sub.concentration1_chr)[0] + "; " + list(meta_sub.agent2)[0] + " " + list(meta_sub.concentration2_chr)[0]
            sns.kdeplot(logdna, ax=ax, label = well + " " + trt, alpha=0.25, bw_adjust=smoothing).set_title(title)
            ax.legend()
        else:
            sns.kdeplot(logdna, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(title)
        plt.close()
    return(fig)

def kde_dna_plot_cell_line(barcode, cell_line, n_wells=None, well_start=0, title = "", 
                       add_legend=False, smoothing=1, column="ldr"):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    wells = get_wells(barcode, cell_line)
    print(wells)
    if n_wells is not None: wells = wells[well_start:(well_start+n_wells)]
    if title == "": title = cell_line + " " + barcode
    fig = kde_dna_plot_wells(barcode, wells, title = title, add_legend=add_legend, 
                         smoothing=smoothing, column = column)
    return(fig)

def get_logdna(barcode, well, remove_na=False, return_dna_also=False, return_dna_colname=True):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    df = read_and_rename_well_data(barcode, well, silent = True)
    dna_orig = df['dna'].copy()
    dna = df['dna'].copy()
    if not remove_na:
        dna[dna < 0] = float('nan')
        logdna = np.log10(dna)
    else:
        dna = dna[dna > 0]
        logdna = np.log10(dna)
        logdna = logdna[ [not x for x in np.isnan(logdna)] ]
        logdna = logdna[ [not x for x in np.isinf(logdna)] ]
    dna_colname = set(df.dna_colname)
    if len(dna_colname) == 1:
        dna_colname = list(dna_colname)[0]
    else:
        print("warning: more than one dna column name: " + barcode + " " + well + " " + str(dna_colname))
        dna_colname = list(dna_colname)
    if not return_dna_colname:
        if return_dna_also:
            return(list(logdna), dna_orig)
        else:
            return(list(logdna))
    else:
        if return_dna_also:
            return(list(logdna), dna_orig, dna_colname)
        else:
            return(list(logdna), dna_colname)


def kde_dna_plot_avg(barcode, cell_line, smoothing=1):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    wells = get_wells(barcode, cell_line)
    #if n_wells is not None: wells = wells[0:n_wells]
    logdna_all = []
    #ax = axs[row,col]
    fig, ax = plt.subplots()
    for well in wells:
        logdna, dna_colname = get_logdna(barcode, well)
        logdna_all.extend(logdna)
    #return(logdna_all)
    sns.kdeplot(logdna_all, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(barcode + " " + cell_line)
    plt.close()
    return(fig)

def kde_dna_plot_all_avg(cell_line, n_barcodes=None, smoothing=1):
    barcodes = get_all_plates_for_cell_line(cell_line)
    fig, ax = plt.subplots()
    if n_barcodes is not None: barcodes = barcodes[0:n_barcodes]
    print(len(barcodes))
    for i in range(len(barcodes)):
        #print(i)
        barcode = barcodes[i]
        if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
        wells = get_wells(barcode, cell_line)
        #if n_wells is not None: wells = wells[0:n_wells]
        logdna_all = []
        for well in wells:
            logdna, dna_colname = get_logdna(barcode, well)
            logdna_all.extend(logdna)
        #return(logdna_all)
        sns.kdeplot(logdna_all, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(cell_line)
        plt.close()
    return(fig)

def kde_dna_plot_plate(barcode, n_wells = None, output_dir="", filename="test_kde_dna.pdf",
                   add_dna_lines=False, add_median_dna_lines=False, smoothing=1):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    cell_lines = get_cell_lines_on_plate(barcode)
    pdf_full = os.path.join(output_dir, filename)
    nb_plots = len(cell_lines)
    if nb_plots != 6: print(barcode + ": " + str(nb_plots) + " cell lines")
    if nb_plots == 0: return(None)
    ncols = 2
    nrows = int(np.ceil(nb_plots/2))
    nrows = 3
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(9, 4*nrows),
                        layout="constrained", sharex= "all")
    #ldr_cutoffs = []
    for row in range(nrows):
        for col in range(ncols):
            i = row*ncols + col
            if i == nb_plots: break
            cell_line = cell_lines[i]
            wells = get_wells(barcode, cell_line)
            if n_wells is not None: wells = wells[0:n_wells]
            logdna_all = []
            ax = axs[row,col]
            dna_colnames = []
            for well in wells:
                logdna, dna_colname = get_logdna(barcode, well)
                dna_colnames.append(dna_colname)
                logdna_all.extend(logdna)
                sns.kdeplot(logdna, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(cell_line)
                #ldr_gates, ldr_lims = dcf_int.get_ldrgates(np.array([10**x for x in logdna]))
                #ldr_cutoff = ldr_gates[1]
                #ldr_cutoffs.append(ldr_cutoff)
                #if add_ldr_line:
                #    ax.axvline(ldr_cutoff,ymin=0, ymax=0.1, color = "red", linestyle = "--")
            sns.kdeplot(logdna_all, ax=ax,  alpha=1, color="black", linewidth=2).set_title(cell_line)
            dna_colname = set(dna_colnames)
            if len(dna_colname) > 1:
                print("warning: more than one dna column name: " + barcode + " " + str(dna_colname))
                xmin = 4
                xmax = 7
            else:
                dna_colname = list(dna_colname)[0]
                if "Hoechst" in dna_colname:
                    xmin = 2.5
                    xmax = 4.5
                else:
                    xmin = 4
                    xmax = 7
            x_lims = (xmin, xmax)
            plt.xlim(x_lims)
            x_ticks = np.arange(np.ceil(xmin), np.floor(xmax)+1)
            plt.xticks(x_ticks)
            ax.tick_params(labelbottom=True)
            #if add_median_ldr_line:
            #    med_cutoff = np.median(ldr_cutoffs)
            #    #print(med_cutoff)
            #    ax.axvline(x=med_cutoff, ymin=0, ymax=1, color = "orange")
    fig.suptitle(barcode)
    plt.savefig(pdf_full)
    return(None)

def kde_dna_plot_all_plates(smoothing=1):
    if not 'folder_dict' in globals(): define_folder_dict('folder_dict')
    for date in folder_dict.keys():
        print(date)
        plates = get_barcodes(date)
        plates.sort()
        for plate in plates:
            folder = "dna_density_plots"
            if not os.path.exists(folder): os.makedirs(folder)
            file = plate + ".pdf"
            if os.path.exists(os.path.join(folder, file)):
                print(plate)
                print("pdf already written")
            else:
                kde_dna_plot_plate(plate, output_dir = folder, filename = file, smoothing=smoothing)
                print(plate)
    return(None)

def kde_dna_replot_plates(smoothing=1):
    hdf = pd.read_csv("hoechst_plates.csv")
    plates = list(hdf['barcode'].copy())
    for plate in plates:
        folder = "dna_density_plots"
        file = plate + ".pdf"
        kde_dna_plot_plate(plate, output_dir = folder, filename = file, smoothing=smoothing)
        print(plate)
    return(None)

def check_edu_and_ph3():
    barcodes = sorted(bc_to_num.values())
    col_dict = {}
    edu = []
    dna = []
    ph3 = []
    hoechst = []
    for barcode in barcodes:
        #print(barcode)
        cell_line = get_cell_lines_on_plate(barcode)[0]
        wells = get_wells(barcode, cell_line)
        for well in wells:
            df = read_well_data(barcode, well)
            if isinstance(df, pd.DataFrame):
                if df.shape[1] > 0:
                    col_dict[barcode] = df.columns
                    edu_list = ["edu" in col.lower() for col in df.columns]
                    edu_ind = np.where(edu_list)
                    if len(edu_ind[0]) > 0:
                        edu_col = df.columns[ edu_ind[0][0] ]
                    else:
                        edu_col = ""
                    ph3_list = ["ph3" in col.lower() for col in df.columns]
                    ph3_ind = np.where(ph3_list)
                    if len(ph3_ind[0]) > 0:
                        ph3_col = df.columns[ ph3_ind[0][0] ]
                    else:
                        ph3_col = ""
                    dna_list = ["dnacontent" in col.lower() for col in df.columns]
                    dna_ind = np.where(dna_list)
                    if len(dna_ind[0]) > 0:
                        dna_col = df.columns[ dna_ind[0][0] ]
                    else:
                        dna_col = ""
                    hoechst_list = ["hoechst" in col.lower() for col in df.columns]
                    hoechst_ind = np.where(hoechst_list)
                    if len(hoechst_ind[0]) > 0:
                        hoechst_col = df.columns[ hoechst_ind[0][0] ]
                    else:
                        hoechst_col = ""
                    edu.append(edu_col)
                    dna.append(dna_col)
                    ph3.append(ph3_col)
                    hoechst.append(hoechst_col)
                    break
                else:
                    pass
            else:
                pass
    df_out = pd.DataFrame({'barcode': barcodes,
         'edu': edu,
         'dna': dna,
         'ph3': ph3,
         'hoechst': hoechst
        })
    return(df_out)


def get_dna_gating_test(dna, ldrint, ldr_gates, x_dna=None, ax=None):
    """Computes gating to claissfy live/dead cells based on DNA content
    
    Parameters
    ----------
    dna : 1d array
         DNA content of cells in a given well
    ldrtxt : 1d array
        ldr txt feature across all cells in a well
    ldr_gates : list of floats
    x_dna : 1d array
        Expected distribution of DNA content (used as x-axis grid)
    ax : subplot object
        provides positional reference for master plot

    Returns
    -------
    dna_gates : list of floats
        inner and outer gates to classify live/dead cells
    """
    logint = np.log10(ldrint)
    logint[np.isnan(logint)] = -10 #dummy value
    
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    log_dna = dcf_int.compute_log_dna(dna, x_dna)
    f_dna = ccg.findpeaks.get_kde(np.array(log_dna), x_dna)

    log_dna_low_ldr = log_dna[ (ldr_gates[1] >= logint) &
                               (logint >= ldr_gates[0])]
    f_dna_low_ldr = ccg.findpeaks.get_kde(log_dna_low_ldr, x_dna)

    g1_loc = dcf_int.get_g1_location(log_dna, x_dna, ldrint, ldr_gates)
    log_dna_g2_range = log_dna[(log_dna > (g1_loc + 0.4 * np.log10(2))) &
                               (ldr_gates[1] >= logint) &
                               (logint >= ldr_gates[0])]

    try:
        f_dna_g2_range = ccg.findpeaks.get_kde(log_dna_g2_range, x_dna)
        g1_g2_pos = dcf_int.get_g1_g2_position(log_dna, x_dna, ldrint, ldr_gates)
        g1_loc = g1_g2_pos[0]
        g2_loc = g1_g2_pos[1]
        dna_gates = [a + b for a, b in zip(
            [g1_g2_pos[i] for i in [0, 0, 1, 1]],
            [(g2_loc-g1_loc) * s for s in [-1.5, -.9, 1.3, 2.2]]
        )]
        y_vals = [np.max(f_dna) * y for y in [0, 1.02, 1.02, 0]]
        inner_x_vals = [dna_gates[i] for i in [1, 1, 2, 2]]
        outer_x_vals = [dna_gates[i] for i in [0, 0, 3, 3]]
        dna_lims = dcf_int.get_dnalims(log_dna, x_dna)
        dna_lims = [np.min((dna_lims[0], dna_gates[0]-0.1)),
                    np.max((dna_lims[1], dna_gates[3]+0.1))]
        return np.array(dna_gates)
    except ValueError:
        return None

def get_g1_location_test(log_dna, x_dna, ldrint, ldr_gates):
    """Computes  ocation of G1 based on DNA content

    Parameters
    ----------
    log_dna : 1d array
       log DNA content of cells in a given well
    x_dna : 1d array
       Expected distribution of DNA content (used as x-axis grid)
    ldrtxt : 1d array
        ldr txt feature across all cells in a well
    ldr_gates : list of floats

    Returns
    -------
    g1_loc : float
       G1 location on log DNA axis
    """
    logint = np.log10(ldrint)
    logint[np.isnan(logint)] = -10 #dummy value
    
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    # Only consider susbet of cells with LDR within ldr_gates
    log_dna_low_ldr = log_dna[(ldr_gates[1] >= logint) &
                              (logint >= ldr_gates[0])]
    f_dna_low_ldr = get_kde(log_dna_low_ldr, x_dna)
    dna_peaks_amp, dna_peaks_loc, _ = findpeaks(f_dna_low_ldr.tolist())
    # Remove lesser peaks
    dna_peaks_loc = dna_peaks_loc[dna_peaks_amp > np.max(dna_peaks_amp/10)]
    dna_peaks_amp = dna_peaks_amp[dna_peaks_amp > np.max(dna_peaks_amp/10)]
    xdna_loc = x_dna[dna_peaks_loc[:4]]  # take the 4 highest peaks
    # compute dna density surrounding peaks
    dna_density = [np.mean(np.array(log_dna > (x - 0.2 * np.log10(2))) &
                           np.array(log_dna < (x + 1.2 * np.log10(2))))
                   for x in xdna_loc] + dna_peaks_amp
    # Find G1 peak
    if len(xdna_loc) == 2:
        g1_loc = np.min(xdna_loc)
    else:
        g1_loc = xdna_loc[np.argmax(dna_density)]
    return g1_loc



## Notes:
# 1) write function to loop through all plates and get all counts -- write each to csv, then to parquet at the end
#### While code is running --- 
# 1) add dna gates to output data frame
# 2) format peak_props as a data frame for output
# 3) write code to write ldr density + gate plots to pdf (later)
# 4) write code to make dna density + gate plots (later)


define_bc_to_num('bc_to_num')
load_well_metadata('meta')
define_regating_df('regate_df')
define_folder_dict('folder_dict')