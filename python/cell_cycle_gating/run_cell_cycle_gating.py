import matplotlib
#matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import logging
from matplotlib.backends.backend_pdf import PdfPages
from cell_cycle_gating import dead_cell_filter as dcf
from cell_cycle_gating import dead_cell_filter_ldrint as dcf_int
from cell_cycle_gating import cellcycle_phases as cc
from cell_cycle_gating import ph3_filter as pf
from cell_cycle_gating.findpeaks import get_kde, findpeaks, get_prominence_reference_level
import seaborn as sns


def run(data, ndict, dfm=None,
        ph3_channel=True, ldr_channel=True,
        px_edu=None, x_ldr=None, control_based_gating=False,
        control_gates=None, fudge_gates=np.array([0, 0, 0, 0]),
        system=None, header=7, ldr_gates=None,
        nuclei_population_name="Nuclei Selected"):
    """Executes cell cycle gating on all wells for which object level
    data is available in object_level_directory. Plots and saves summary pdf
    of DNA v EdU distribution with automated gatings. A dataframe summarizing
    results for each well is also saved.

    Parameters
    ----------
    object_level_directory : str
                           name of folder containing object level data for a single plate
    dfm : Optional[pandas dataframe]
        metadata table of experimental design. Default is None.
    ph3_channel : Optional[bool]
                True if data on pH3 intensity is in object level data.
    nuclei_population_name : Optional[str]
        Name of the nuclei "population" in your Harmony analysis protocol. Used
        to identify which data files to read.

    Returns
    -------
    df : pandas dataframe
       well level summary of number of live/dead cells and fraction of
       cells in each phase of the cell cycle
    """
    if not os.path.isdir('results'):
        os.mkdir('results')
    if not os.path.isdir('logs'):
        os.mkdir('logs')
    plt.ioff()
    df_summary = pd.DataFrame()
    identity_dict = {}
    df_gates = pd.DataFrame()
    
    if os.path.isdir(data):
        logfile = "logs/%s.log" % data.split('[')[0]        
        if dfm is not None:
            dfm_ord = merge_metadata(dfm, data, nuclei_population_name)
            if control_based_gating:
                dfm_ord = dfm_ord[dfm_ord.agent == 'DMSO'].copy()
            object_level_data = dfm_ord['object_level_file'].tolist()
        else:
            object_level_data = [s for s in os.listdir(data)
                                 if f'{nuclei_population_name}[0].txt' in s]
            dfm_ord=None
    else:
        logfile = "logs/%s.log" % data.split('.txt')[0]
        df_input = pd.read_table(data, header=header)
        df_input = df_input.rename(columns=ndict)
        if dfm is not None:
            barcode = data.split('.txt')[0]
            dfm = dfm[dfm.barcode == barcode].copy()
            metadata_wells = dfm.well.unique()
            df_input['well'] = df_input['well'].astype("category")
            #df_input.well.cat.set_categories(metadata_wells, inplace=True)
            df_input['well'] = df_input.well.cat.set_categories(metadata_wells)
            df_input = df_input.sort_values(['well'])
            dfm_ord = dfm.sort_values(['cell_line', 'agent', 'concentration'])
            dfm_ord.index = dfm_ord['well']
            if control_based_gating:
                dfm_ord = dfm_ord[dfm_ord.agent == 'DMSO'].copy()
            object_level_data = dfm_ord.well.unique()    
        else:
            dfm_ord = None
            object_level_data = df_input.well.unique()

    logger = logging.getLogger()
    if (logger.hasHandlers()):
        logger.handlers.clear()
    logger.setLevel(logging.DEBUG)
        
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    # tell the handler to use this format
    console.setFormatter(formatter)
    logger.addHandler(console)

    errorfile = logging.FileHandler(logfile)
    errorfile.setLevel(logging.ERROR)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    errorfile.setFormatter(formatter)
    logger.addHandler(errorfile)
    #logging.basicConfig(filename=logfile, level=logging.ERROR)

                            
    nb_plots = len(object_level_data)    
    if control_based_gating:
        pdf_pages = PdfPages('results/control_summary_%s.pdf' % data.split('.txt')[0])
    else:
        pdf_pages = PdfPages('results/summary_%s.pdf' % data.split('.txt')[0])
        
    nb_plots_per_page = 10
    # nb_pages = int(np.ceil(nb_plots / float(nb_plots_per_page)))
    for i, file in enumerate(object_level_data):
        if i % nb_plots_per_page == 0:
            fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
      
        try:
            if os.path.isdir(data):
                df = pd.read_table('%s/%s' % (data, file))  
                well = re.search('result.(.*?)\[', file).group(1)
                well = "%s%s" % (well[0], well[1:].zfill(2))
                df['well'] = well
            else:
                df = df_input[df_input.well == file].copy()
                df['ldrint'] = df['Cell: LDRrawINT (DDD-bckgrnd)'] - \
                    df['Cell: LDRbackground (DDD-bckgrnd)']
                well = file
                if 'dna' not in df.columns.tolist():
                    df['dna'] = df['Cell: Average Intensity_Average (DDD)'].multiply(
                        df['Cell: Area_Average (DDD)'])
                if ph3_channel and 'ph3' not in df.columns.tolist():
                    df['ph3'] = df['Cell: pH3rawINT (DDD-bckgrnd)'] - \
                        df['Cell: pH3background (DDD-bckgrnd)']    
    
            df = map_channel_names(df, ndict)
            # if system=='ixm':
            #     edu = np.array(df['edu'].tolist())
            #     if ph3_channel:
            #         ph3 = np.array(df['ph3'].tolist())
            #         cells_notnan = ~(np.isnan(edu) | np.isnan(ph3))
            #     else:
            #         cells_notnan = ~np.isnan(edu)
            #     ldr = np.array(df['ldr'].tolist())
            #     ldr = ldr[cells_notnan]
            #     ldr_min = np.max((100, ldr.min() - 100))
            #     ldr_max = ldr.max() + 100
            #     x_ldr = np.arange(ldr_min, ldr_max, 100)
            fractions, gates, cell_identity = gate_well(df, dfm_ord=dfm_ord,
                                                        ph3_channel=ph3_channel,
                                                        ldr_channel=ldr_channel,
                                                        px_edu=px_edu, x_ldr=x_ldr,
                                                        control_based_gating=control_based_gating,
                                                        control_gates=control_gates,
                                                        fudge_gates=fudge_gates,
                                                        fig=fig, plot_num=i,
                                                        ldr_gates=ldr_gates)
            #df_summary = df_summary.append(fractions, ignore_index=True)
            #df_gates = df_gates.append(gates, ignore_index=True)
            df_summary = pd.concat([df_summary, pd.DataFrame([fractions])], ignore_index=True)
            df_gates = pd.concat([df_gates, pd.DataFrame([gates])], ignore_index=True)
            identity_dict[well] = cell_identity
       
        except ValueError as err:
            logging.error("%s in well %s" % (err, well))
            #pass
        except TypeError as err:
            logging.error("%s in well %s" % (err, well))
            #pass
        except IndexError as err:
            logging.error("%s in well %s" % (err, well))
            #pass
        # your code that will (maybe) throw
        except np.linalg.LinAlgError as err:
            if 'Singular matrix' in str(err):
                logging.error("%s in well %s" % (err, well))
            #pass
        except ZeroDivisionError as err:
            logging.error("%s in well %s" % (err, well))
            #pass
        #except pd.io.common.EmptyDataError as err:
        except pd.errors.EmptyDataError as err:
            logging.error("%s in well %s" % (err, well))
            #pass
        if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
            plt.tight_layout()
            pdf_pages.savefig(fig)
            logging.info("Completed analysis for %d out of %d wells" %
                         (i+1, len(object_level_data)))
            plt.close('all')
    pdf_pages.close()
    summary_cols = ['well', 'cell_count__total',
                    'mean_Sphase_edu',
                    # 'cell_count', 'cell_count__dead',
                    'G1', 'S', 'G2',
                    'S_dropout', 'subG1', 'beyondG2']
    if ph3_channel:
        summary_cols.append('M')
    if ldr_channel:
        summary_cols += ['cell_count', 'cell_count__dead']
    df_summary = df_summary[summary_cols]
    # Merge summary table with metadata if provided
    if dfm is not None:
        df_summary.index = df_summary['well'].tolist()
        df_summary = pd.concat([dfm_ord, df_summary], axis=1)
        df_summary = df_summary.loc[:, ~df_summary.columns.duplicated()]
    # Merge summary table with corpse count if available
    if os.path.isdir(data):
        dfc = get_corpse_count(data)
    else:
        dfc = None
    if dfc is not None:
        df_summary.index = df_summary['well'].tolist()
        df_summary = pd.concat([df_summary, dfc], axis=1)
    if control_based_gating:
        df_summary = df_summary[df_summary.agent == 'DMSO'].copy()
        df_summary.to_csv('results/control_summary_%s.csv' % data.split('.txt')[0], index=False)
        return df_summary, df_gates
    else:
        df_summary.to_csv('results/summary_%s.csv' % data.split('.txt')[0], index=False)        
        return df_summary  #, identity_dict
    

def gate_well(df, dfm_ord=None, ph3_channel=True, ldr_channel=True,
              px_edu=None, x_ldr=None, control_based_gating=False,
              control_gates=None, fudge_gates=np.array([0, 0, 0, 0]),
              fig=None, plot_num=0, ldr_gates=None, ldr_control_cutoff=2):
    """Gating on a single well
    
    Parameters
    ----------
    df : pandas.DataFrame
       single cell data for a single well of the plate

    Return
    ------
    fractions : dict
       Fraction of cells classified to be in each phase of the cell cycle
    gates : dict
       Value of EdU and DNA at gates for given cell line and well
    cell_identity : list[int]
       mapping of each cell to a phase of the cell cycle
    """
    well = df.well.unique()[0]
    well = "%s%s" % (well[0], well[1:].zfill(2))
    edu = np.array(df['edu'].tolist())
    dna = np.array(df['dna'].tolist())
    #edu = edu.astype(np.float)
    edu = edu.astype(float)
    #edu[edu < 0] = np.nan

    if ph3_channel:
        ph3 = np.array(df['ph3'].tolist())
        cells_notnan = ~(np.isnan(edu) | np.isnan(ph3))
        ph3 = ph3[cells_notnan]
    else:
        cells_notnan = ~np.isnan(edu)
    edu = edu[cells_notnan]
    #if edu.min() < 0:
    #    edu += np.abs(edu.min())
    dna = dna[cells_notnan]
         
    # Get phases based on DNA and EdU
    if dfm_ord is not None:
        # dfm_ord.index = dfm_ord.well
        cell_line = dfm_ord.loc[well, 'cell_line']
        agent = dfm_ord.loc[well, 'agent']
        conc = dfm_ord.loc[well, 'concentration']
        title = "%s %s %s (%s um)" % (well, cell_line, agent, conc)
    else:
        title = well
    if control_gates is not None:
       control_dna_gates = control_gates[control_gates.cell_line == cell_line][
           'dna_gates'].mean()
       #control_edu_gates = control_gates[control_gates.cell_line == cell_line][
       #    'edu_gates'].mean()
       control_edu_gates=None
       control_ldr_cutoff =  control_gates[control_gates.cell_line == cell_line]['ldr_cutoff'].mean()
       control_edu_cutoff_bounds = [control_gates[control_gates.cell_line == cell_line]['mean_Gphase_edu'].mean(),
                                   control_gates[control_gates.cell_line == cell_line]['mean_Sphase_edu'].mean()]
                                   

    else:
        control_dna_gates = None
        control_edu_gates = None
        control_ldr_cutoff = 1.2
        control_edu_cutoff_bounds = None
        #control_edu_gates = rb_gates

    # Get live dead
    if ldr_channel:
        ldrint = np.array(df['ldrint'].tolist())
        ldrint = ldrint[cells_notnan]
        logldr_peak = dcf_int.get_ldr_peak_val(ldrint)
        ldr_gates, ldr_lims = dcf_int.get_ldrgates(ldrint, ldr_control_cutoff, peak_loc=control_ldr_cutoff)
        logint = np.log10(ldrint)
        logint[np.isnan(logint)] = -10 
        ldr_inner = ((ldr_gates[1] >= logint) & (logint >= ldr_gates[0]))
        if np.sum(ldr_inner) < 50:
            #dna = None
            dna_gates=None
        else:
            #dna = df['dna']
            dna_gates = dcf_int.get_dna_gating(dna, ldrint, ldr_gates)
            # if dna_gates is None:
            #     dna=None
        cell_fate_dict, outcome = dcf_int.live_dead(ldrint, ldr_gates=ldr_gates, dna=dna,
                                                    dna_gates= dna_gates,
                                                    ldr_control_cutoff=ldr_control_cutoff)
        outcome = np.array(outcome)
        live_cols = [s for s in list(cell_fate_dict.keys()) if 'alive' in s]
        dead_cols = [s for s in list(cell_fate_dict.keys()) if 'dead' in s]
        a = 0
        d = 0
        for col in live_cols:
            a += cell_fate_dict[col]
        for col in dead_cols:
            d += cell_fate_dict[col]
    else:
        outcome = np.array([1] * len(dna))
    ### if all cells are dead
    if sum(outcome>=1) == 0:
        fractions = {}
        gates = {}
        cell_identity = outcome
        cols = []
        if ph3_channel:
            cols = ['subG1', 'G1', 'S', 'S_dropout', 'G2', 'M', 'beyondG2']
        else:
            cols = ['subG1', 'G1', 'S', 'S_dropout', 'G2', 'beyondG2']
        for x in cols:
            fractions[x] = 0
        fractions['mean_Sphase_edu'] = float('nan')
    else:
        try:
            fractions, cell_identity, gates = cc.plot_summary(dna[outcome>=1], edu[outcome>=1], fig=None,
                                                              title=title,
                                                              plot='scatter',
                                                              plot_num=plot_num,
                                                              px_edu=px_edu,
                                                              control_dna_gates=control_dna_gates,
                                                              control_edu_gates=control_edu_gates,
                                                              fudge_gates=np.array(fudge_gates),
                                                              control_edu_cutoff_bounds=control_edu_cutoff_bounds)
        except ValueError:
            fractions = {}
            gates = {}
            cell_identity = np.array([0]*len(dna[outcome>=1]))
            
    
        # Revise phases based on pH3
        if ph3_channel:
            #grid_size = (5, 2)
            #plot_num = plot_num % 10
            #rel_pos = np.unravel_index(plot_num, grid_size)
            #axp = plt.subplot2grid(grid_size, rel_pos)
            
            log_ph3 = pf.compute_log_ph3(ph3[outcome>=1])
            #x_ph3 = np.arange(2.5, 8, 0.02)
            #f_ph3 = get_kde(log_ph3, x_ph3)
            #peak_amp, peak_loc, peak_width = findpeaks(f_ph3.tolist(), npeaks=2, thresh=0.1)
            #if len(peak_amp) >= 2:
            #    cutoff_loc, _ = get_prominence_reference_level(
            #        f_ph3.tolist(), peak_amp[1], peak_loc[1])
            #    ph3_cutoff = x_ph3[cutoff_loc]
            #del_ph3 = x_ph3[peak_loc[0]] - log_ph3.min()
            #ph3_cutoff = x_ph3[peak_loc[0]] + del_ph3
            #else:
            f_ph3, ph3_cutoff, ph3_lims = pf.get_ph3_gates(ph3[outcome>=1], cell_identity)
            fractions, cell_identity = pf.evaluate_Mphase(log_ph3, ph3_cutoff, cell_identity)
            #sns.kdeplot(log_ph3, ax=axp)
            #axp.vlines(ph3_cutoff, 0, 2)
            #axp.text(ph3_cutoff+0.1, 1, str(fractions['M']))
            #axp.set_title(title, fontsize=6)
    
        edu_live = edu[outcome >=1]   
        fractions['mean_Sphase_edu'] = edu_live[cell_identity==2].mean()
        gates['mean_Sphase_edu'] = np.log10(edu_live[cell_identity==2].mean())
        gates['mean_Gphase_edu'] = np.log10(edu_live[cell_identity==1].mean())
    
    if ldr_channel:
        fractions['cell_count'] = a
        fractions['cell_count__dead'] = d
        gates['ldr_cutoff'] = logldr_peak
    fractions['well'] = well
    fractions['cell_count__total'] = len(dna)

    if dfm_ord is not None:
        gates['well'] = well
        gates['cell_line'] = cell_line


    return fractions, gates, cell_identity


def merge_metadata(dfm, obj, nuclei_population_name):
    """Merges metadata with object level file and
    sorts files by cell line, drug and concetration.

    Parameters
    ----------
    dfm : pandas dataframe
        metadata file with wells mapped to treatment conditions
    obj : str
        path to folder containing object level data
    nuclei_population_name : str
        Name of the nuclei "population" in the Harmony analysis protocol.

    Returns
    -------
    dfmc : pandas dataframe
        mapping between wells, treatment, and filename
        containing object level data.
        Sorted by cell line, drug and concentration
    """
    # Merge agent and concentration columns in the case of
    # combination experiments.
    dfm = process_metadata_file(dfm)

    barcode = obj.split('[')[0]
    dfm = dfm[dfm.barcode == barcode].copy()
    dfm.index = dfm.well.tolist()
    wells = []
    object_level_files = [s for s in os.listdir(obj)
                          if f'{nuclei_population_name}[0].txt' in s]
    for file in object_level_files:
        well = re.search('result.(.*?)\[', file).group(1)
        well = "%s%s" % (well[0], well[1:].zfill(2))
        wells.append(well)
    dfmap = pd.DataFrame(list(zip(wells, object_level_files)),
                         columns=['well', 'object_level_file'])
    dfmap.index = dfmap.well.tolist()
    dfmc = pd.concat([dfm, dfmap], axis=1)
    dfmc = dfmc.loc[:, ~dfmc.columns.duplicated()]
    dfmc = dfmc.dropna(subset=['object_level_file'])
    dfmc = dfmc.sort_values(['cell_line', 'agent', 'concentration'])
    return dfmc


def process_metadata_file(dfmeta):
    """Merge agent and concentration columns in the case of
       combination experiments.

    Parameters
    ----------
    dfmeta : pandas dataframe
       metadata file with wells mapped to treatment condition.

    Returns
    -------
    dfm : pandas dataframe
       metadata file with agent and concentration columns merged in the case of
    combination experiments.
    """
    dfm = dfmeta.copy()
    agent_cols = [a for a in dfm.columns.tolist() if 'agent' in a]
    if len(agent_cols) > 1:
        dfm[agent_cols] = dfm[agent_cols].replace(['', ' '], np.nan)
        dfm['agent'] = dfm[agent_cols].apply(
            lambda x: ','.join(x[x.notnull()]), axis=1)
        conc_cols = [a for a in dfm.columns.tolist() if 'concentration' in a]
        dfm[conc_cols] = dfm[conc_cols].astype(str)
        dfm[conc_cols] = dfm[conc_cols].replace(['0.0', '0'], np.nan)
        dfm['concentration'] = dfm[conc_cols].apply(
            lambda x: ','.join(x[x.notnull()]), axis=1)
    return dfm


def get_corpse_count(obj):
    """Return corpse count dataframe

    Parameters
    ----------
    obj : str
        path to folder containing object level data


    Returns
    -------
    dfc : pandas dataframe
        datatable mapping well to corpse counts
    """
    corpse_file = [f for f in os.listdir(obj)
                   if 'CorpseCountOnly' in f and '.txt' in f]
    if corpse_file:
        df_corpse = pd.read_table("%s/%s" %
                                  (obj, corpse_file[0]))
        wells = ["%s%s" % (w[0], w[1:].zfill(2))
                 for w in df_corpse.WellName.tolist()]
        df_corpse.index = wells
        nd = {'Corpses - Number of Objects': 'corpse_count'}
        df_corpse = df_corpse.rename(columns=nd)
        df_corpse = df_corpse['corpse_count'].copy()
        dfc = pd.DataFrame(df_corpse)
        return dfc
    else:
        return None


def map_channel_names(df, ndict):
    if 'well' in df.columns.tolist():
        ndict = {key:val for key, val in ndict.items() if val != 'well'}
    df = df.rename(columns=ndict)
    df = df.loc[:, ~df.columns.duplicated()]
    return df


def get_ixmc_barcode(filename, corpse_col=None):
    """
    Parameters
    ----------
    filename : str
      Name of .txt file containing well-level DDD data
    
    Returns
    -------
    barcode : str
      Plate barcode 
    """
    headers = []
    f = open(filename)
    for _ in range(8):
        headers.append(f.readline())
    barcode = headers[4].split('=')[1].split('"\n')[0]
    # if '_' not in barcode:
    #     date = barcode[:6]
    #     plate_number = re.search(r'I*(.*)', barcode).group(1)
    #     filler = re.search(r'%s([^"]*)%s' % (date, plate_number)).group(1)
    columns = headers[-1].split('"\t"')
    columns = [s.split('"\n')[0] for s in columns]
    if corpse_col in columns:
        barcode = 'corpse_%s' % barcode
    f.close()    
    return barcode


