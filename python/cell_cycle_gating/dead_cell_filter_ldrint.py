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
    logint = logint[ [not x for x in np.isnan(logint)] ]
    logint = logint[ [not x for x in np.isinf(logint)] ]
    fig, ax = plt.subplots()
    x, y = sns.kdeplot(logint, ax=ax).get_lines()[0].get_data()
    plt.close()
    peak_locs, _ = find_peaks(-y)
    #print(peak_locs)
    cc = x[peak_locs]
    #print(cc)
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



def get_counts(batch, filename, ndict, ldr_control_cutoff=2, is_csv=False, peak_loc = 1.2):
    well = re.search('result.(.*?)\[', filename).group(1)
    well = "%s%s" % (well[0], well[1:].zfill(2))
    if is_csv:
        df = pd.read_csv("%s/%s" % (batch, filename))
    else:
        df = pd.read_table("%s/%s" % (batch, filename))
    barcode = batch.split('[')[0]
    df = df.rename(columns=ndict)
    #ldrint = df['Nuclei Selected - LDRINT']
    #ldrint = df['ldrint']
    ldrint = df['ldr']
    ldr_gates, ldr_lims = get_ldrgates(ldrint, ldr_control_cutoff, peak_loc)
    print(ldr_gates)
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
    logint = logint[ [not x for x in np.isnan(logint)] ]
    logint = logint[ [not x for x in np.isinf(logint)] ]
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

    
