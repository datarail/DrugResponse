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


### functions to read data and call joint-gating, etc.

###### write new function based on Coralie's data here #######

#### plotting functions

def plot_ldr_all_plates(df, barcode_col= 'barcode', cell_line_col = 'cell_line', 
                    well_col = 'well', ldr_col = 'ldr', output_pdf = True,
                   n_wells = None, output_dir="", smoothing=1,
                  add_ldr_line = False, add_median_ldr_line = True,
                  new_gating_algorithm = True, peak_loc = float('-inf')):
    barcodes = df[barcode_col].unique().tolist()
    for barcode in barcodes:
        qq = "barcode == " + "'" + barcode + "'"
        df_tmp = df.query(qq)
        fname = barcode + "_ldr_plot" + ".pdf"
        plot_ldr_plate(df_tmp, barcode, barcode_col=barcode_col, cell_line_col=cell_line_col,
                       well_col=well_col, ldr_col=ldr_col, output_pdf=output_pdf,
                       n_wells=n_wells, output_dir=output_dir, filename=fname,
                       smoothing=smoothing, add_ldr_line=add_ldr_line, 
                       add_median_ldr_line=add_median_ldr_line,
                       new_gating_algorithm=new_gating_algorithm,
                       peak_loc=peak_loc
                       )
    return(None)

def plot_ldr_all_plates_no_gating(df, barcode_col= 'barcode', cell_line_col = 'cell_line', 
                    well_col = 'well', ldr_col = 'ldr', output_pdf = False,
                   n_wells = None, output_dir="", smoothing=1,
                  new_gating_algorithm = True):
    barcodes = df[barcode_col].unique().tolist()
    for barcode in barcodes:
        qq = "barcode == " + "'" + barcode + "'"
        df_tmp = df.query(qq)
        fname = barcode + "_ldr_plot" + ".pdf"
        plot_ldr_plate(df_tmp, barcode, barcode_col=barcode_col, cell_line_col=cell_line_col,
                       well_col=well_col, ldr_col=ldr_col, output_pdf=output_pdf,
                       n_wells=n_wells, output_dir=output_dir, filename=fname,
                       smoothing=smoothing, add_ldr_line=False, 
                       add_median_ldr_line=False,
                       new_gating_algorithm=new_gating_algorithm
                       )
    return(None)

def plot_ldr_plate_no_gating(df, barcode, barcode_col= 'barcode', cell_line_col = 'cell_line', 
                    well_col = 'well', ldr_col = 'ldr', output_pdf = False,
                    n_wells = None, output_dir="", filename="test_ldr_plate.pdf",
                    new_gating_algorithm=True):
    
    plot_ldr_plate(df=df, barcode=barcode, cell_line_col=cell_line_col,
                   well_col=well_col, ldr_col=ldr_col, output_pdf=output_pdf,
                   n_wells=n_wells, output_dir=output_dir, filename=filename,
                   add_ldr_line=False, add_median_ldr_line=False,
                   new_gating_algorithm=new_gating_algorithm)
    return(None)
    

def plot_ldr_plate(df, barcode, barcode_col= 'barcode', cell_line_col = 'cell_line', 
                    well_col = 'well', ldr_col = 'ldr', output_pdf = True,
                    n_wells = None, output_dir="", filename="test_ldr_plate.pdf",
                    add_ldr_line = False, add_median_ldr_line = True,
                    new_gating_algorithm = True,
                    ### options for old gating algorithm (get_ldrgates function)
                    peak_loc = float('-inf'),
                    ### options for new gating algorithm (get_ldrgates_new function)
                    smoothing=1,
                    first_peak_min=float('-inf'),
                    min_prominence=0,
                    min_peak_height=0.02,
                    min_peak_distance=0.5,
                    single_peak_cutoff=float('inf'),
                    mixture_backup_method="valley", silent=True,
                    return_peaks_only=False,
                    ### options for joint gating
                    window_size = 0.2
                  ):
    cell_lines = df[cell_line_col].unique().tolist()
    pdf_full = os.path.join(output_dir, filename)
    nb_plots = len(cell_lines)
    if nb_plots != 6: print(barcode + ": " + str(nb_plots) + " cell lines")
    if nb_plots == 0: return(None)
    ncols = 2
    nrows = int(np.ceil(nb_plots/2))
    nrows = 3
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(9, 4*nrows),
                        layout="constrained", sharex= "all")
    ldrint_all = df[ldr_col].copy()
    ldrint_all = ldrint_all[ldrint_all > 0]
    logint_all = np.log10(ldrint_all)
    logint_all = logint_all.dropna()
    xmin = min(logint_all)
    xmax = max(logint_all)
    for row in range(nrows):
        for col in range(ncols):
            ldr_cutoffs = []
            i = row*ncols + col
            if i == nb_plots: break
            if i >= len(cell_lines): break
            cell_line = cell_lines[i]
            qq = cell_line_col + " == '" + cell_line + "'" + " & " + barcode_col + " == '" + barcode + "'" 
            df_small = df.query(qq)
            if df_small.shape[0] == 0: break
            wells = df_small[well_col].unique().tolist()
            if n_wells is not None: wells = wells[0:n_wells]
            ax = axs[row,col]
            print(cell_line)
            ### initialize a list and dict for outputs from ldr gating function
            df_list = []
            props_dict = {}
            for well in wells:
                qq2 = well_col + " == " + "'" + well + "'"
                df_well = df_small.query(qq2)
                if df_well.shape[0] == 0: break
                ldrint = df_well[ldr_col].copy()
                ldrint = ldrint[ldrint > 0]
                logint = np.log10(ldrint.copy())
                logint = logint.dropna()
                sns.kdeplot(logint, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(cell_line)
                if add_ldr_line or add_median_ldr_line:
                    if not new_gating_algorithm:
                        ldr_gates, ldr_lims = dcf_int.get_ldrgates(np.array([10**x for x in logint]),
                                                                   peak_loc = peak_loc)
                        ldr_cutoff = ldr_gates[1]
                        ldr_cutoffs.append(ldr_cutoff)
                    else:
                        pdict = get_ldrgates_new(ldrint, 
                                        smoothing=smoothing,
                                        show=True,
                                        first_peak_min=first_peak_min,
                                        min_prominence=min_prominence, 
                                        min_peak_height=min_peak_height, min_peak_distance=min_peak_distance,
                                        single_peak_cutoff=single_peak_cutoff,
                                        mixture_backup_method=mixture_backup_method, 
                                        silent=silent,
                                        return_peaks_only=return_peaks_only,
                                        suppress_fig=True)
                        ldr_cutoff = pdict['ldr_cutoff']
                        ldr_cutoffs.append(ldr_cutoff)
                        ### construct results data frame
                        pdict['well'] = well
                        pdict['barcode']=barcode
                        pdict['cell_line']=cell_line
                        cols = ['barcode', 'cell_line', 'ldr_cutoff', 'peak1', 'peak2', 
                                    'peak1_height', 'peak2_height','shelf','method_used', 'ldr_cutoff_mixture', 
                                'ldr_cutoff_valley', 'ldr_cutoff_middle']
                        # cols = ['well','ldr_cutoff', 'peak1', 'peak2', 
                        #             'peak1_height', 'peak2_height','shelf','method_used', 'ldr_cutoff_mixture', 
                        #         'ldr_cutoff_valley', 'ldr_cutoff_middle']
                        dd = {k:v for (k,v) in pdict.items() if k in cols}
                        pdict2 = {k:v for (k,v) in pdict.items() if not k in cols}
                        res_df = pd.DataFrame(data=dd, index=[0])
                        res_df = res_df[cols]
                        df_list.append(res_df)
                        props_dict[well] = pdict2
                if add_ldr_line:
                    ax.axvline(ldr_cutoff,ymin=0, ymax=1, color = "red", linestyle = "--", alpha = 0.1)
            #sns.kdeplot(logint_all, ax=ax, color = "red").set_title(cell_line)
            
            #xmin = -2
            #xmax = 6
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
    if output_pdf:
        print(pdf_full)
        plt.savefig(pdf_full)
    ### organize output dataframe and dictionary
    res_df_full = pd.concat(df_list)
    res_df_full = res_df_full.reset_index(drop=True)
    ldr_cutoff_median = np.nanmedian(res_df_full.ldr_cutoff)
    ldr_cutoff_min = ldr_cutoff_median - window_size
    ldr_cutoff_max = ldr_cutoff_median + window_size
    res_df_full['ldr_cutoff_median'] = ldr_cutoff_median
    res_df_full['ldr_cutoff_min'] = ldr_cutoff_min
    res_df_full['ldr_cutoff_max'] = ldr_cutoff_max
    res_df_full['ldr_cutoff_final'] = np.clip(res_df_full['ldr_cutoff'], ldr_cutoff_min, ldr_cutoff_max)
    #results = {'df': res_df_full, 'peak_props': props_dict}
    return(res_df_full)


def plot_ldr_barcode_cell_line(df, barcode, cell_line, barcode_col= 'barcode', cell_line_col = 'cell_line', 
                               well_col = 'well', ldr_col = 'ldr', add_legend = False, smoothing=1):
    qq = cell_line_col + " == '" + cell_line + "'" + " & " + barcode_col + " == '" + barcode + "'" 
    df_small = df.query(qq)
    wells = df_small[well_col].unique().tolist()
    fig = plot_ldr_wells(df_small, wells, ldr_col = ldr_col, well_col = well_col, add_legend=add_legend, smoothing = smoothing)
    return(fig)

def plot_ldr_wells(df, wells, title = "", add_legend=False, smoothing=1, ldr_col = "ldr", well_col = "well"):
    fig, ax = plt.subplots()
    #ax.set_title(title, fontsize=12)
    for well in wells:
        qq = well_col + " == " + "'" + well + "'"
        df_tmp = df.query(qq)
        ldrint = df_tmp[ldr_col].copy()
        ldrint[ldrint < 0] = float('nan')
        logint = np.log10(ldrint)
        #x, y = sns.kdeplot(logint, ax=ax, color = "grey", alpha=0.5).get_lines()[0].get_data()
        if add_legend:
            sns.kdeplot(logint, ax=ax, label = well, alpha=0.25, bw_adjust=smoothing).set_title(title)
            ax.legend()
        else:
            sns.kdeplot(logint, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(title)
        plt.close()
    return(fig)

### main function (joint gating)
def get_ldrgates_joint(ldrint_list, wells, window_size = 0.2, smoothing = 1, barcode = "", cell_line = ""):
    """Joint gating of many wells, based on ldr intensities
    This function gates each well separately and then 

    Parameters
    ----------
    ldrint_list : list of lists
        Each item in the list should be an ldrint feature (a list itself) across all cells in a well
    wells: list
        a list of the names of each well: e.g. "M07", "C03", etc. -- each name must be unique
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
        pdict = get_ldrgates_new(ldrint.copy(), show=False, silent=True, smoothing=smoothing)
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
        #ldrint = ldrint[ldrint > 0]
        #logint = np.log10(ldrint.copy())
        #logint = logint.dropna()
        #fig, ax = plt.subplots()
        #sns.kdeplot(logint, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(barcode + cell_line)
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

### main function (single-well gating)
def get_ldrgates_new(ldrint, smoothing=1, show=False, first_peak_min=float('-inf'),
                     min_prominence=0, min_peak_height=0.02, min_peak_distance=0.5,
                     single_peak_cutoff=float('inf'),
                     mixture_backup_method="valley", silent=True,
                    return_peaks_only=False, suppress_fig=False):
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
        if not suppress_fig:
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

### main peak finding function
# inputs: x,y from sns.kdeplot output
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


### LDR cutoff finding functions

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


#### helpers:

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

# input x and y from the sns.kdeplot function
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


### extra (not needed, but could be useful for debugging, testing, etc.)

def plot_density_derivative(barcode, well):
    if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
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