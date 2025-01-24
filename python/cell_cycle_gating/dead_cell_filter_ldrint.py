from scipy.signal import find_peaks
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
import os
import pandas as pd
from cell_cycle_gating.findpeaks import get_kde, findpeaks
from cell_cycle_gating import smooth
from scipy.stats.mstats import mquantiles as quantile
import matplotlib.pyplot as plt

from pomegranate.gmm import GeneralMixtureModel
from pomegranate.distributions import *
import torch



def get_ldrgates(ldrint, ldr_control_cutoff=2, peak_loc=1.2):
    """Gating based on ldr intensities

    Parameters
    ----------
    ldrint : 1d array
        ldr txt feature across all cells in a well
    ldr_cutoff : float
        default cutoff if gating fails.

    Returns
    -------
    ldr_gates : list of floats
    ldr_lims : list of floats
        limits of ldr inensity feature that defines x_lims for plots
    """
    ldrint = ldrint[ldrint > 0]
    logint = np.log10(ldrint)
    logint = logint[~np.isnan(logint)]
    logint = logint[~np.isinf(logint)]
    fig, ax = plt.subplots()
    x, y = sns.kdeplot(logint, ax=ax).get_lines()[0].get_data()
    plt.close()
    peak_locs, _ = find_peaks(-y)
    cc = x[peak_locs]
    try:
        ldr_cutoff = cc[cc > peak_loc][0]
    except IndexError:
        ldr_cutoff = np.quantile(logint, 0.99)
        #ldr_cutoff = ldr_control_cutoff
    ldr_gates = np.array([-np.inf, ldr_cutoff])
    ldr_lims = np.array([x.min(), x.max()])
    return ldr_gates, ldr_lims


def compute_log_dna(dna, x_dna=None):
    """Computes log of DNA content bounded by x_dna[2], x_dna[-3]

    Parameters
    ----------
    dna : 1D array
        DNA content of cells in a given well
    x_dna : 1D array
        Expected distribution of DNA content (used as x-axis grid)

    Return
    ------
    log_dna : 1D array
        log transformed DNA content
    """
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    dna_upper_bound = 10 ** x_dna[-3]
    dna_lower_bound = 10 ** x_dna[2]
    dna_upper_bounded = [d if d < dna_upper_bound else dna_upper_bound
                         for d in dna]
    dna_bounded = [d if d > dna_lower_bound else dna_lower_bound
                   for d in dna_upper_bounded]
    log_dna = np.array([np.log10(d) for d in dna_bounded])
    return log_dna



def get_g1_location(log_dna, x_dna, ldrint, ldr_gates):
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


def get_g2_location(log_dna, x_dna, ldrint, ldr_gates, g1_loc):
    """Computes location of G2 based on DNA content

    Parameters
    ----------
    log_dna : 1d array
       log DNA content of cells in a given well
    x_dna : 1d array
       Expected distribution of DNA content (used as x-axis grid)
    ldrtxt : 1d array
        ldr txt feature across all cells in a well
    ldr_gates : list of floats
    g1_loc : numpy float
        G1 location on log DNA scale

    Returns
    -------
    g2_loc : numpy float
        G2 location on log DNA scale
    """
    # Get G2 peak and location
    # Only consider subset of cells witt LDR internsity within ldr_gates and
    # DNA content > (g1_loc + 0.4 * log10(2))
    logint = np.log10(ldrint)
    logint[np.isnan(logint)] = -10 #dummy value

    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    log_dna_g2_range = log_dna[(log_dna > (g1_loc + 0.4 * np.log10(2))) &
                               (ldr_gates[1] >= logint) &
                               (logint >= ldr_gates[0])]
    f_dna_g2_range = get_kde(log_dna_g2_range, x_dna)
    f_smooth = smooth.smooth(f_dna_g2_range, 5, 'flat')
    peak_amp, peak_loc, _ = findpeaks(f_smooth.tolist())
    peak_loc = peak_loc[peak_amp > np.max(peak_amp/10)]
    xdna_loc = x_dna[peak_loc]
    xdna_loc = xdna_loc[xdna_loc > (g1_loc + 0.5 * np.log10(2))]
    if len(xdna_loc) > 1:
        g2_loc = xdna_loc[np.argmin(
            np.abs((xdna_loc - (g1_loc + np.log10(2))))
        )]
    elif len(xdna_loc) == 1:
        g2_loc = xdna_loc[0]
    else:
        g2_loc = g1_loc + np.log10(2)
    return g2_loc


def get_g1_g2_position(log_dna, x_dna, ldrint, ldr_gates):
    """Wrapper function that returns G1 and G2 location
    based on log DNA content

    Parameters
    ----------
    log_dna : 1d array
       log DNA content of cells in a given well
    x_dna : 1d array
       Expected distribution of DNA content (used as x-axis grid)
    ldrint : 1d array
        ldr intensiy feature across all cells in a well
    ldr_gates : list of floats

    Returns
    -------
    g1_g2_pos : list of floats
        G1 and G2 location on log DNA scale
    """
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    g1_loc = get_g1_location(log_dna, x_dna, ldrint, ldr_gates)
    g2_loc = get_g2_location(log_dna, x_dna, ldrint, ldr_gates, g1_loc)
    g1_g2_pos = [g1_loc, g2_loc]
    return g1_g2_pos


def get_dnalims(log_dna, x_dna=None):
    """ Outer bounds on DNA content to use as x_lim for plots

    Parameters
    ----------
    log_dna : 1d array
        log DNA content of cells in a given well
    x_dna : 1d array
        Expected distribution of DNA content (used as x-axis grid)

    Returns
    -------
    dna_lims : list of floats
    """
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    dna_lims = (quantile(log_dna, [5e-3, 0.995]) +
                [(2.5 * (x_dna[1] - x_dna[0])) * x for x in [-1, 1]])
    return dna_lims


def get_dna_gating(dna, ldrint, ldr_gates, x_dna=None, ax=None):
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
    log_dna = compute_log_dna(dna, x_dna)
    f_dna = get_kde(np.array(log_dna), x_dna)

    log_dna_low_ldr = log_dna[ (ldr_gates[1] >= logint) &
                               (logint >= ldr_gates[0])]
    f_dna_low_ldr = get_kde(log_dna_low_ldr, x_dna)

    g1_loc = get_g1_location(log_dna, x_dna, ldrint, ldr_gates)
    log_dna_g2_range = log_dna[(log_dna > (g1_loc + 0.4 * np.log10(2))) &
                               (ldr_gates[1] >= logint) &
                               (logint >= ldr_gates[0])]

    try:
        f_dna_g2_range = get_kde(log_dna_g2_range, x_dna)
        g1_g2_pos = get_g1_g2_position(log_dna, x_dna, ldrint, ldr_gates)
        g1_loc = g1_g2_pos[0]
        g2_loc = g1_g2_pos[1]
        dna_gates = [a + b for a, b in zip(
            [g1_g2_pos[i] for i in [0, 0, 1, 1]],
            [(g2_loc-g1_loc) * s for s in [-1.5, -.9, 1.3, 2.2]]
        )]
        y_vals = [np.max(f_dna) * y for y in [0, 1.02, 1.02, 0]]
        inner_x_vals = [dna_gates[i] for i in [1, 1, 2, 2]]
        outer_x_vals = [dna_gates[i] for i in [0, 0, 3, 3]]
        dna_lims = get_dnalims(log_dna, x_dna)
        dna_lims = [np.min((dna_lims[0], dna_gates[0]-0.1)),
                    np.max((dna_lims[1], dna_gates[3]+0.1))]
        return np.array(dna_gates)
    except ValueError:
        return None



def live_dead(ldrint, ldr_gates=None,
              dna=None, dna_gates=None,
              x_dna=None, ax=None, ldr_control_cutoff=2):
    """Assign classification to individual cells as live/dead based on
    ldrint and DNA content.
    If ax is not None, plots pie chart of fraction live/dead
    1. alive = selected+others, where selected is within
                 inner DNA gate and within LDR
    2. dead = anything outside of DNA outer gating and LDR gating
    3. total = alive + dead; selected + others + dead

    Parameters
    ----------
    ldrint : 1d array
       ldr int feature across all cells in a well
    ldr_gates : list of floats
    dna : 1d array
       DNA content of cells in a given well
    dna_gates : list of floats
    ldr_cutoff : float
       default cutoff if
    x_dna : 1d array
       Expected distribution of DNA content (used as x-axis grid)
    ax : subplot object

    Returns
    -------
    alive : int
       number of cells classied as alive
    dead : int
       numer of cells classified as dead
    outcome : 1d array
       classification of each cell as live(>=0) or dead (-1).
       should have same length as ldrtxt
    """
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)

    outcome = [0] * len(ldrint)
    logint = np.log10(ldrint)
    logint[np.isnan(logint)] = -10 #dummy value
    
    if ldr_gates is  None:
        ldr_gates, _ = get_ldrgates(ldrint, ldr_control_cutoff)
    ldr_outer = (logint < ldr_gates[0]) | (logint > ldr_gates[1])
    outcome = [-1 if b else 0 for b in ldr_outer]
    #dead = np.sum([1 for ot in outcome if ot == -1])
    alive = np.sum([1 for ot in outcome if ot >= 0])
    selected = 'DNA information unavailable'
    others = 'DNA information unavailable'
    dead_ldrpos = np.sum(ldr_outer)
    cell_fate_dict = {'alive': alive, 'alive_subg1': 0, 'alive_beyondg2': 0,
                          'dead_ldrpos': dead_ldrpos, 'dead_subg1': 0}

    if dna_gates is not None:
        log_dna = compute_log_dna(dna, x_dna)
        dna_outermost = (log_dna < dna_gates[0]) | (log_dna > dna_gates[3])
        dead_ldrpos = np.sum(ldr_outer)
        dead_subg1 = np.sum((ldr_outer==False) & (log_dna < dna_gates[0]))
        alive_beyondg2 = np.sum((ldr_outer==False) & (log_dna > dna_gates[2]))
        alive_subg1 = np.sum((ldr_outer==False) & (log_dna > dna_gates[0]) & (log_dna < dna_gates[1]))
        dna_inner = ((log_dna > dna_gates[1]) &
                     (log_dna < dna_gates[2]) &
                     (ldr_outer==False))
        alive = np.sum(dna_inner)
        #outcome = [-1 if d else 1 if s else 0
        #           for d, s in zip((ldr_outer | dna_outermost), dna_inner)]
        outcome = ((1 * dna_inner) # normal live cells
                   + (1.5 * ((ldr_outer==False) & (log_dna > dna_gates[2]))) # live but higher than G2
                   + (-1 * ((ldr_outer==False) & (log_dna < dna_gates[0]))) # dead very low G1
                   + (1.25 * ((ldr_outer==False) & (log_dna > dna_gates[0]) & (log_dna < dna_gates[1]))) # alive lower than G1
                   + (-2 * ldr_outer)) 
        cell_fate_dict = {'alive': alive, 'alive_subg1': alive_subg1, 'alive_beyondg2': alive_beyondg2,
                          'dead_ldrpos': dead_ldrpos, 'dead_subg1': dead_subg1}
        #alive = np.sum([1 for ot in outcome if ot >= 0])
        #dead = np.sum([1 for s in outcome if s == -1])
        #selected = np.sum([1 for s in outcome if s == 1])
        #others = np.sum([1 for s in outcome if s == 0])
        if ax is not None:
            ax.pie([alive, alive_subg1, alive_beyondg2, dead_ldrpos, dead_subg1],
                   labels=['alive', 'alive_subg1', 'alive_beyondg2', 'dead_ldrpos', 'dead_subg1'],
                   explode=(0.1, 0.1, 0.1, 0.1, 0.1), autopct='%1.1f%%')
            ax.axis('equal')
    else:
        if ax is not None:
            ax.pie([alive, dead_ldrpos], labels=['alive', 'dead_ldrpos'],
                   explode=(0.1, 0.1), autopct='%1.1f%%')
            ax.axis('equal')
    return cell_fate_dict, outcome



def get_counts(batch, filename, ndict, ldr_control_cutoff=2):
    well = re.search('result.(.*?)\[', filename).group(1)
    well = "%s%s" % (well[0], well[1:].zfill(2))
    df = pd.read_table("%s/%s" % (batch, filename))
    barcode = batch.split('[')[0]
    df = df.rename(columns=ndict)
    #ldrint = df['Nuclei Selected - LDRINT']
    ldrint = df['ldrint']
    ldr_gates, ldr_lims = get_ldrgates(ldrint, ldr_control_cutoff)
    logint = np.log10(ldrint)
    logint[np.isnan(logint)] = -10 
    ldr_inner = ((ldr_gates[1] >= logint) & (logint >= ldr_gates[0]))
    if np.sum(ldr_inner) < 50:
        dna = None
        dna_gates=None
    else:
        dna = df['dna']
        dna_gates = get_dna_gating(dna, ldrint, ldr_gates)
        if dna_gates is None:
            dna=None
    cell_fate_dict, outcome = live_dead(ldrint, ldr_gates=ldr_gates, dna=dna, dna_gates= dna_gates,
                                        ldr_control_cutoff=ldr_control_cutoff)
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
    dfs['cell_count'] = a
    dfs['cell_count__dead'] = d
    dfs['cell_count__total'] = len(ldrint)
    return(dfs)


def batch_run(batch, ndict):
    files = [s for s in os.listdir(batch) if s.endswith('Nuclei Selected[0].txt')]
    dfl = []
    for f in files:
        dfl.append(get_counts(batch, f, ndict))
    dfc = pd.concat(dfl)
    dfc.index = range(len(dfc))
    return(dfc)

               
def get_ldr_peak_val(ldrint):
    #ldrint = df['ldrint']
    ldrint = ldrint[ldrint > 0]
    logint = np.log10(ldrint)
    logint = logint[~np.isnan(logint)]
    logint = logint[~np.isinf(logint)]
    fig, ax = plt.subplots()
    x, y = sns.kdeplot(logint, ax=ax).get_lines()[0].get_data()
    plt.close()
    peak_locs, _ = find_peaks(y)
    max_loc = np.argmax(y[peak_locs])
    peak_val = x[peak_locs][max_loc]
    return(peak_val)
                       

def summary_peak_val(batch, ndict):
    files = [s for s in os.listdir(batch) if s.endswith('Nuclei Selected[0].txt')]
    dfl = []
    for f in files:
        dfl.append(get_ldr_peak_val(batch, f, ndict))
    dfc = pd.concat(dfl)
    dfc.index = range(len(dfc))
    return(dfc)

###### Functions for new gating methods below (written by Nick Clark, January 2025) ----------------

def get_counts_plate_batch(df,
                           barcode_col_design= 'barcode', 
                           cell_line_col_design = 'cell_line', 
                           well_col_design = 'well',
                           ldr_col_data = 'ldrint', 
                           well_col_data = 'Well Name',
                           dna_col_data = 'Cell: DNAcontent (DD-bckgrnd)',
                           output_pdf = True,
                   n_wells = None, output_dir="", smoothing=1,
                  add_ldr_line = False, add_median_ldr_line = True,
                  new_gating_algorithm = True, metadata_cols = [],
                  peak_loc = float('-inf'), main_dir = "", well_file_ending = 'test].csv'):
    """Get live/dead cell counts for all plates and wells listed in a design layout file.
    The function loops over the plates and gates them one by one. It considers all wells of the same cell line on the same plate to be a group, first finding gates for them individually, then selecting the median of these gates as the final gate for the group.

    Parameters
    ----------
    df : pd.DataFrame
       This is a dataframe with the design layout (i.e. metadata) for the plates that need to be gated. It must have columns listing the plate barcodes (e.g. "201117_combo_12"), well names (e.g. "A01"), and cell_line (e.g. "HCC70").
    barcode_col_design : str, optional
        The name of the column in the design dataframe (df) with the plate barcodes. Default is "barcode".
    cell_line_col_design : str, optional
        The name of the column in the design dataframe (df) with the cell line names. Default is "cell_line".
    well_col_design : str, optional
        The name of the column in the design dataframe (df) with the well names. Default is "well".
    ldr_col_data : str, optional
        The name of the column in the well-level data with Live/Dead Red fluorescence values. This is passed on to the "get_counts_plate" function. Default is "ldrint".
    well_col_data: str, optional
        The name of the column in the well-level data with well names. This is passed on to the "get_counts_plate" function. Default is "Well Name".
    dna_col_data: str, optional
        The name of the column in the well-level data with DNA fluorescence values. This is passed on to the "get_counts_plate" function. Default is "Cell: DNAcontent (DD-bckgrnd)".
    output_pdf : bool, optional
        Whether to output a pdf with plots of gating for each well or not. Default is True.
    n_wells : int, optional
        If None, then all wells will be gated, otherwise only the first "n_wells" will be gated. This is included mostly for testing purposes. Default is None (gating all wells).
    output_dir : str, optional
        Directory for the output. Default is the current directory.
    smoothing : int, optional
        Amount of "smoothing" to use for the kernel density function when gating. Default is 1.
    add_ldr_line : bool, optional
        Whether to include lines to show Live/Dead Red gating for each individual well. This usually looks messy, so default is False.
    add_median_ldr_line : bool, optional
        Whether to include a line showing the median Live/Dead Red gating across wells in a group (same cell line, same plate). Default is True.
    new_gating_algorithm : bool, optional
        Whether to use the updated gating algorithm (implemented by Nick Clark, Jan. 2025) or not. This algorithm makes the gating of individual wells much more robust for edge cases (e.g. small blip in the density function, "shelf" rather than "valley" between live and dead peaks) and finds an optimal gate by fitting a Gaussian mixture model. To make the gating even more robust, it then finds the median gate for each group of wells (same cell line, same plate) and uses this gate for the group. Default is True.
    metadata_cols : list of str, optional
        Other metadata columns in the design dataframe. Default is an empty list: [].
    peak_loc : float
        Fluorescence value to start at when looking for peaks. Setting this parameter can sometimes help when gating functions find "false peaks" or blips in the density function before the true live/dead peaks. Default is float('-inf'), i.e. this parameter is not used by default. This parameter was more useful with previous gating methods, but with the updated gating methods it is not necessary.
    main_dir : str
        Directory with LDR Fluorescence data for each plate. The structure should be [main_dir]/[barcode]/[well_fluorescence_csv_files]. Default is the current directory.
    well_file_ending : str
        File ending that identifies files with data for each well. Default is 'test].csv'.
        
    Returns
    -------
        A pandas dataframe with the gates and live/dead cell counts for each well for all plates in the batch.
    """
    barcodes = df[barcode_col_design].unique().tolist()
    df_list = []
    for barcode in barcodes:
        qq = "barcode == " + "'" + barcode + "'"
        df_tmp = df.query(qq)
        df = df[df[cell_line_col_design].notna()]
        fname = barcode + "_ldr_plot" + ".pdf"
        plate_dir = os.path.join(main_dir, barcode)
        df_tmp = get_counts_plate(df_tmp, barcode,
            barcode_col_design=barcode_col_design,
            cell_line_col_design=cell_line_col_design,
                       well_col_design=well_col_design, ldr_col_data=ldr_col_data,
                       well_col_data = well_col_data,
                       dna_col_data = dna_col_data,
                       output_pdf=output_pdf,
                       n_wells=n_wells, output_dir=output_dir, filename=fname,
                       smoothing=smoothing, add_ldr_line=add_ldr_line, 
                       add_median_ldr_line=add_median_ldr_line,
                       new_gating_algorithm=new_gating_algorithm,
                       metadata_cols = metadata_cols,
                       peak_loc=peak_loc,
                       plate_data_dir = plate_dir,
                       well_file_ending = well_file_ending
                       )
        df_list.append(df_tmp)
    df_full = pd.concat(df_list)
    df_full = df_full.reset_index(drop=True)
    return(df_full)

def get_counts_plate(df, barcode, plate_data_dir=None, 
                    well_file_ending = 'test].csv',
                    barcode_col_design= 'barcode', 
                    cell_line_col_design = 'cell_line', 
                    well_col_design = 'well', 
                    ldr_col_data = 'ldrint',
                    well_col_data = 'Well Name',
                    dna_col_data = 'Cell: DNAcontent (DD-bckgrnd)',
                    output_pdf = True,
                    n_wells = None, output_dir="",
                    add_ldr_line = False, add_median_ldr_line = True,
                    new_gating_algorithm = True,
                    metadata_cols = [],
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
                    window_size = 0,
                    testing = False
                  ):
    """Get live/dead cell counts for all wells on a single plate.
    The function considers all wells of the same cell line on the same plate to be a group, first finding gates for them individually, then selecting the median of these gates as the final gate for the group.

    Parameters
    ----------
    df : pd.DataFrame
       This is a dataframe with the design layout (i.e. metadata) for the plate to be gated. It must have columns listing the plate barcodes (e.g. "201117_combo_12"), well names (e.g. "A01"), and cell_line (e.g. "HCC70").
    barcode : str
        The plate barcode.
    plate_data_dir : str, optional
        The directory with csv files of LDR fluorescence data for individual wells. Default is None, which means the function will look in [currect directory]/[barcode] for the data files, selecting any that are named ending in 'test].csv'.
    well_file_ending : str
        File ending that identifies files with data for each well. Default is 'test].csv'.
    barcode_col_design : str, optional
        The name of the dataframe column with the plate barcodes. Default is "barcode".
    cell_line_col_design : str, optional
        The name of the dataframe column with the cell line names. Default is "cell_line".
    well_col_design : str, optional
        The name of the dataframe column with the well names. Default is "well".
    ldr_col_data : str, optional
        The name of the column in the well-level data with Live/Dead Red fluorescence values. Default is "ldrint".
    well_col_data: str, optional
        The name of the column in the well-level data with well names. Default is "Well Name".
    dna_col_data: str, optional
        The name of the column in the well-level data with DNA fluorescence values. Default is "Cell: DNAcontent (DD-bckgrnd)".
    output_pdf : bool, optional
        Whether to output a pdf with plots of gating for each well or not. Default is True.
    n_wells : int, optional
        If None, then all wells will be gated, otherwise only the first "n_wells" will be gated. This is included mostly for testing purposes. Default is None (gating all wells).
    output_dir : str, optional
        Directory for the output. Default is the current directory.
    smoothing : int, optional
        Amount of "smoothing" to use for the kernel density function when gating. Default is 1.
    add_ldr_line : bool, optional
        Whether to include lines to show Live/Dead Red gating for each individual well. This usually looks messy, so default is False.
    add_median_ldr_line : bool, optional
        Whether to include a line showing the median Live/Dead Red gating across wells in a group (same cell line, same plate). Default is True.
    new_gating_algorithm : bool, optional
        Whether to use the updated gating algorithm (implemented by Nick Clark, Jan. 2025) or not. This algorithm makes the gating of individual wells much more robust for edge cases (e.g. small blip in the density function, "shelf" rather than "valley" between live and dead peaks) and finds an optimal gate by fitting a Gaussian mixture model. To make the gating even more robust, it then finds the median gate for each group of wells (same cell line, same plate) and uses this gate for the group. Default is True.
    metadata_cols : list of str, optional
        Other metadata columns in the design dataframe. Default is an empty list: [].
    peak_loc : float
        Fluorescence value to start at when looking for peaks. Setting this parameter can sometimes help when gating functions find "false peaks" or blips in the density function before the true live/dead peaks. Default is float('-inf'), i.e. this parameter is not used by default. This parameter was more useful with previous gating methods, but with the updated gating methods it is not necessary.
        
    Returns
    -------
        A pandas dataframe with the gates and live/dead cell counts for each well for all plates in the batch. Prints a pdf with gating for each group (cell line) on the plate. Note: It is written write now to assume up to 6 cell lines per plate... may error if there are more.
    """
    if plate_data_dir is None:
        plate_data_dir = barcode
    filename = barcode + "_LDR_gating.pdf"
    qq = "barcode == " + "'" + barcode + "'"
    df = df.query(qq)
    df = df[df[cell_line_col_design].notna()]
    cell_lines = df[cell_line_col_design].unique().tolist()
    pdf_full = os.path.join(output_dir, filename)
    nb_plots = len(cell_lines)
    if nb_plots != 6: print(barcode + ": " + str(nb_plots) + " cell lines")
    if nb_plots == 0: return(None)
    ncols = 2
    nrows = int(np.ceil(nb_plots/2))
    nrows = 3
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(9, 4*nrows),
                        layout="constrained", sharex= "all")
    # wells = df[well_col_data].unique().tolist()
    # for well in wells:
    #     df_tmp = read_well_data(barcode, well, plate_data_dir,
    #                                      ldr_col_data = ldr_col_data,
    #                                      well_col_data = well_col_data,
    #                                      dna_col_data = dna_col_data)
    # ldrint_all = df['ldr'].copy()
    # ldrint_all = ldrint_all[ldrint_all > 0]
    # logint_all = np.log10(ldrint_all)
    # logint_all = logint_all.dropna()
    # xmin = min(logint_all)
    # xmax = max(logint_all)
    df_list_full = []
    xmin = float('inf')
    xmax = float('-inf')
    for row in range(nrows):
        for col in range(ncols):
            ldr_cutoffs = []
            i = row*ncols + col
            if i == nb_plots: break
            if i >= len(cell_lines): break
            cell_line = cell_lines[i]
            qq = cell_line_col_design + " == '" + cell_line + "'" + " & " + barcode_col_design + " == '" + barcode + "'" 
            df_small = df.query(qq)
            if df_small.shape[0] == 0: break
            wells = df_small['well'].unique().tolist()
            if n_wells is not None: wells = wells[0:n_wells]
            ax = axs[row,col]
            print(cell_line)
            ### initialize a list and dict for outputs from ldr gating function
            df_list = []
            #props_dict = {}
            #well_files = [s for s in os.listdir(plate_data_dir) if s.endswith(well_file_ending)]
            ### read into a data frame with
            #df_plate_data = read_plate_data(barcode, well_file_ending)
            for well in wells:
                #qq2 = well_col + " == " + "'" + well + "'"
                #df_well = df_small.query(qq2)
                #well_file = [f for f in well_files if well in f]
                #df_well = pd.read_csv(well_file[0])
                df_well = read_well_data(barcode, well, plate_data_dir,
                                         ldr_col_data = ldr_col_data,
                                         well_col_data = well_col_data,
                                         dna_col_data = dna_col_data)
                if df_well.shape[0] == 0: break
                ldrint = df_well['ldr'].copy()
                ldrint = ldrint[ldrint > 0]
                logint = np.log10(ldrint.copy())
                logint = logint.dropna()
                xmin_tmp = min(logint)
                xmax_tmp = max(logint)
                xmin = min(xmin_tmp, xmin)
                xmax = max(xmax_tmp, xmax)
                sns.kdeplot(logint, ax=ax, alpha=0.25, bw_adjust=smoothing).set_title(cell_line)
                if add_ldr_line or add_median_ldr_line:
                    if not new_gating_algorithm:
                        ldr_gates, ldr_lims = get_ldrgates(np.array([10**x for x in logint]),
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
                        pdict['barcode'] = barcode
                        pdict['cell_line'] = cell_line
                        cols = ['well','barcode', 'cell_line', 'ldr_cutoff']
                        # cols = ['well','ldr_cutoff', 'peak1', 'peak2', 
                        #             'peak1_height', 'peak2_height','shelf','method_used', 'ldr_cutoff_mixture', 
                        #         'ldr_cutoff_valley', 'ldr_cutoff_middle']
                        dd = {k:v for (k,v) in pdict.items() if k in cols}
                        for col in metadata_cols:
                            len1 = df_well[col].unique().size
                            if len1 == 1:
                                dd[col] = df_well[col].unique()[0]
                            elif len1 < 1:
                                dd[col] = None
                                msg = "Warning for barcode: " + barcode + " well: " + well
                                msg = msg + "no values for metadata column " + "'" + col + "'"
                                warnings.warn(msg)
                            elif len1 > 1:
                                dd[col] = None
                                msg = "Warning for barcode: " + barcode + " well: " + well
                                msg = msg + "multiple value for metadata column " + "'" + col + "'"
                                msg = msg + ' '.join(str(x) for x in df_well[col].unique())
                                warnings.warn(msg)
                        #pdict2 = {k:v for (k,v) in pdict.items() if not k in cols}
                        ### dataframe with ldr cutoff (one row)
                        tmp_df = pd.DataFrame(data=dd, index=[0])

                        #### get counts from ldr gates
                        ldr_gates = np.array([-np.inf, ldr_cutoff])
                        dna_gates = None
                        dna = None
                        cell_fate_dict, outcome = live_dead(ldrint, ldr_gates=ldr_gates, dna=dna, dna_gates= dna_gates)
                        live_cols = [s for s in list(cell_fate_dict.keys()) if 'alive' in s]
                        dead_cols = [s for s in list(cell_fate_dict.keys()) if 'dead' in s]
                        a = 0
                        d = 0
                        for col in live_cols:
                            a += cell_fate_dict[col]
                        for col in dead_cols:
                            d += cell_fate_dict[col]
                        dfs = pd.DataFrame(cell_fate_dict, index=[0])
                        #print('dna_gates')
                        #print(dna_gates)
                        dfs['cell_count'] = a
                        dfs['cell_count__dead'] = d
                        dfs['cell_count__total'] = len(ldrint)
                        tmp_df_full = pd.concat([tmp_df, dfs], axis = 1)
                        df_list.append(tmp_df_full)
                        #props_dict[well] = pdict2
                if add_ldr_line:
                    ax.axvline(ldr_cutoff,ymin=0, ymax=1, color = "red", linestyle = "--", alpha = 0.1)
            ### organize output dataframe and dictionary ##
            res_df = pd.concat(df_list)
            res_df = res_df.reset_index(drop=True)
            ldr_cutoff_median = np.nanmedian(res_df.ldr_cutoff)
            ldr_cutoff_min = ldr_cutoff_median - window_size
            ldr_cutoff_max = ldr_cutoff_median + window_size
            res_df['ldr_cutoff_median'] = ldr_cutoff_median
            res_df['ldr_cutoff_min'] = ldr_cutoff_min
            res_df['ldr_cutoff_max'] = ldr_cutoff_max
            res_df['ldr_cutoff_final'] = np.clip(res_df['ldr_cutoff'], ldr_cutoff_min, ldr_cutoff_max)
            df_list_full.append(res_df)
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
    res_df_full = pd.concat(df_list_full)
    res_df_full = res_df_full.reset_index(drop=True)
    if new_gating_algorithm:
        func_settings = ['smoothing', 'first_peak_min', 'min_prominence', 'min_peak_height',
                     'min_peak_distance', 'single_peak_cutoff', 'mixture_backup_method', 'peak_loc']
        res_df_full['gating_algorithm'] = "new algorithm"
    else:
        func_settings = ['peak_loc']
        res_df_full['gating_algorithm'] = "old algorithm"
    if not testing:
        cell_count_cols = ['cell_count', 'cell_count__dead', 'cell_count__total']
        phase_cols = ['alive', 'alive_subg1', 'alive_beyondg2', 'dead_ldrpos', 'dead_subg1']
        cols = ['well', 'barcode', 'cell_line'] + metadata_cols + ['ldr_cutoff_final'] + phase_cols + cell_count_cols
        res_df_full = res_df_full.filter(items = cols)
        res_df_full.rename(columns={'ldr_cutoff_final':'ldr_cutoff'}, inplace=True)
    
    #results = {'df': res_df_full, 'peak_props': props_dict}
    return(res_df_full)

### main function (single-well gating)
def get_ldrgates_new(ldrint, smoothing=1, show=False, first_peak_min=float('-inf'),
                     min_prominence=0, min_peak_height=0.02, min_peak_distance=0.5,
                     single_peak_cutoff=float('inf'),
                     mixture_backup_method="valley", silent=True,
                    return_peaks_only=False, suppress_fig=False):
    ldrint = ldrint.copy()
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

###### helper functions -----------

## Get the filename for well-level intensities for a given barcode and well
## input:
##    barcode: a plate barcode, e.g. '210406_combo_71'
##    well: a well of interest, e.g. 'D06'
## output:
##    full path/filename of the well-level data
def get_well_file(barcode, well, data_dir = ''):
    #date = date_from_barcode(barcode)
    #data_dir = get_data_dir(barcode)
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
def read_well_data(barcode, well,
                   plate_data_dir = None,
                   ldr_col_data = 'ldrint',
                   well_col_data = 'Well Name',
                   dna_col_data = 'Cell: DNAcontent (DD-bckgrnd)'):
    if plate_data_dir is None:
        plate_data_dir = barcode
    ff = get_well_file(barcode, well, plate_data_dir)
    col_dict = {}
    col_dict[ldr_col_data] = 'ldr'
    col_dict[well_col_data] = 'well'
    col_dict[dna_col_data] = 'dna'
    df = pd.read_csv(ff)
    df = df.rename(columns=col_dict)
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
    #if len(barcode) in [2,3]: barcode = barcode_from_number(barcode)
    df = read_well_data(barcode, well)
    df = rename_df_columns(df, silent = silent, hoechst_as_dna=hoechst_as_dna)
    return(df)