from cell_cycle_gating.findpeaks import get_kde, findpeaks
import numpy as np
from scipy.stats.mstats import mquantiles as quantile
import matplotlib.pyplot as plt
# from itertools import compress
from cell_cycle_gating import smooth
from scipy.stats import gaussian_kde
import matplotlib.gridspec as gridspec
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42



def get_ldrgates(ldrtxt, x_ldr=None):
    """Gating based on ldr intensities

    Parameters
    ----------
    ldrtxt : 1d array
        ldr txt feature across all cells in a well
    x_ldr : 1d array
        uniformly distributed 1d grid based on expected
        range of ldr txt

    Returns
    -------
    ldr_gates : list of floats
       gating based on minima of kernel density estimate of ldr txt
    """
    if x_ldr is None:
        mx = np.max(ldrtxt.tolist())+0.01
        x_ldr = np.arange(-0.01, mx, 0.0002)
    f_ldr = get_kde(ldrtxt, x_ldr)  # ldrtxt should be an array
    peak_amp, peak_loc, peak_width = findpeaks(f_ldr.tolist(), npeaks=1)

    # Find location of minimun on the right
    f_neg = [-x for x in f_ldr[peak_loc[0]:]]
    _, trough_loc, _ = findpeaks(f_neg, npeaks=1)
    # If peakfinder cannot find peak minima, use ldrwidth_5x as default
    if np.any(trough_loc):
        trough_loc = trough_loc[0] + peak_loc[0] - 1
    else:
        trough_loc = peak_loc + (5 * peak_width[0])

    # choose LDR cutoff based on half-proximal width and right trough of peak
    ldrwidth_5x = peak_loc + (5 * peak_width[0])
    ldrwdith_2p5 = peak_loc + (2.5 * peak_width[0])
    cutoff_index_1 = len(x_ldr) - 2
    cutoff_index_2 = np.max([3,
                             np.min([trough_loc, ldrwidth_5x]),
                             ldrwdith_2p5])

    ldr_cutoff = x_ldr[np.min([cutoff_index_1, int(cutoff_index_2)])]

    ldr_gates = [-np.inf, ldr_cutoff]
    return np.array(ldr_gates)


def get_ldrlims(ldrtxt, x_ldr=None):
    """Limits of ldr txt feature that define x_lims for plots

    Parameters
    ----------
    ldrtxt : 1d array
        ldr txt feature across all cells in a well
    x_ldr : 1d array
        uniformly distributed 1d grid based on expected
        range of ldr txt

    Returns
    -------
    ldr_lims : list of floats
        limits of ldr txt feature that define x_lims for plots
    """
    if x_ldr is None:
        mx = np.max(ldrtxt.tolist())+0.01
        x_ldr = np.arange(-0.01, mx, 0.0002)
    ldr_lims = (quantile(ldrtxt, [5e-3, 0.995]) +
                [(2.5 * (x_ldr[1] - x_ldr[0])) * x for x in [-1, 1]])
    return ldr_lims


def plot_ldr_gating(ldrtxt, x_ldr=None, ldr_gates=None,
                    ldr_lims=None, ax=None):
    """Summary plot of gating based on gating based on LDR intensities

    Parameters
    ----------
    ldrtxt : 1d array
       ldr txt feature across all cells in a well
    x_ldr : 1d array
       uniformly distributed 1d grid based on expected
       range of ldr txt
    ldr_gates : list of floats
       min and max gating based on ldr txt
    ldr_lims : list of floats
       outer bouns of ldr txt to set as x_lim for pltos
    ax : plot obj
       provides positional reference for master plot

    Returns
    -------
    """
    if x_ldr is None:
        mx = np.max(ldrtxt.tolist())+0.01
        x_ldr = np.arange(-0.01, mx, 0.0002)
    f_ldr = get_kde(ldrtxt, x_ldr)
    if not ldr_gates:
        ldr_gates = get_ldrgates(ldrtxt, x_ldr)
    if not ldr_lims:
        ldr_lims = get_ldrlims(ldrtxt, x_ldr)
    log_frac = np.log10(f_ldr+np.max(f_ldr)/100) - np.log10(np.max(f_ldr)/100)

    if ax is None:
        ax = plt.figure()
    ax.plot(x_ldr, log_frac)
    x_vals = [ldr_gates[1],
              np.max([ldr_gates[0], np.min(x_ldr)]),
              np.max([ldr_gates[0], np.min(x_ldr)]),
              ldr_gates[1], ldr_gates[1]]
    y_vals = [np.log10(np.max(f_ldr)) * y for y in [0, 0, 0.5, 0.5, 0]]
    ax.plot(x_vals, y_vals, 'r')
    ax.set_xlim(ldr_lims)
    f_ldr_max = np.log10(np.max(f_ldr)) - np.log10(np.max(f_ldr)/100) + 0.1
    ax.set_ylim([0, f_ldr_max])
    ax.set_xlabel('LDRtxt intensity')
    ax.set_ylabel('kernel density estimate')
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


def get_g1_location(log_dna, x_dna, ldrtxt, ldr_gates):
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
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    # Only consider susbet of cells with LDR within ldr_gates
    log_dna_low_ldr = log_dna[(ldr_gates[1] >= ldrtxt) &
                              (ldrtxt >= ldr_gates[0])]
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


def get_g2_location(log_dna, x_dna, ldrtxt, ldr_gates, g1_loc):
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
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    log_dna_g2_range = log_dna[(log_dna > (g1_loc + 0.4 * np.log10(2))) &
                               (ldr_gates[1] >= ldrtxt) &
                               (ldrtxt >= ldr_gates[0])]
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


def get_g1_g2_position(log_dna, x_dna, ldrtxt, ldr_gates):
    """Wrapper function that returns G1 and G2 location
    based on log DNA content

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
    g1_g2_pos : list of floats
        G1 and G2 location on log DNA scale
    """
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    g1_loc = get_g1_location(log_dna, x_dna, ldrtxt, ldr_gates)
    g2_loc = get_g2_location(log_dna, x_dna, ldrtxt, ldr_gates, g1_loc)
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


def get_dna_gating(dna, ldrtxt, ldr_gates, x_dna=None, ax=None):
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
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    log_dna = compute_log_dna(dna, x_dna)
    f_dna = get_kde(np.array(log_dna), x_dna)

    log_dna_low_ldr = log_dna[(ldr_gates[1] >= ldrtxt) &
                              (ldrtxt >= ldr_gates[0])]
    f_dna_low_ldr = get_kde(log_dna_low_ldr, x_dna)

    g1_loc = get_g1_location(log_dna, x_dna, ldrtxt, ldr_gates)
    log_dna_g2_range = log_dna[(log_dna > (g1_loc + 0.4 * np.log10(2))) &
                               (ldr_gates[1] >= ldrtxt) &
                               (ldrtxt >= ldr_gates[0])]
    f_dna_g2_range = get_kde(log_dna_g2_range, x_dna)

    g1_g2_pos = get_g1_g2_position(log_dna, x_dna, ldrtxt, ldr_gates)
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
    if ax is not None:
        ax.plot(x_dna, f_dna_low_ldr, '--r')
        ax.plot(x_dna, f_dna, '-k')
        ax.plot(x_dna, f_dna_g2_range, ':')
        ax.plot(inner_x_vals, y_vals, '-r', linewidth=2)
        ax.plot(outer_x_vals, y_vals, '-r')
        ax.set_xlabel('log10 (DNA content)')
        ax.set_ylabel('kernel density estimate')
        ax.set_xlim(dna_lims)
    return np.array(dna_gates)


def plot_ldr_dna_scatter(dna, ldrtxt, x_dna=None, x_ldr=None, ax=None):
    """Plot of LDR and DNA scatter with associated gates
    
    Parameters
    ----------
    dna : 1d array
        DNA content of cells in a given well
    ldrtxt : 1d array
        ldr txt feature across all cells in a well
    x_dna : 1d array
        Expected distribution of DNA content (used as x-axis grid)
    x_ldr : 1d array
        uniformly distributed 1d grid based on expected
        range of ldr txt
    ax : subplot object
        provides positional reference for master plot
    
    Returns
    -------
    """
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    if x_ldr is None:
        mx = np.max(ldrtxt.tolist())+0.01
        x_ldr = np.arange(-0.01, mx, 0.0002)
    log_dna = compute_log_dna(dna, x_dna)
    xy = np.vstack([log_dna, ldrtxt])
    z = gaussian_kde(xy)(xy)

    ldr_gates = get_ldrgates(ldrtxt, x_ldr)

    g1_g2_pos = get_g1_g2_position(log_dna, x_dna, ldrtxt, ldr_gates)
    g1_loc = g1_g2_pos[0]
    g2_loc = g1_g2_pos[1]
    dna_gates = [a + b for a, b in zip(
        [g1_g2_pos[i] for i in [0, 0, 1, 1]],
        [(g2_loc-g1_loc) * s for s in [-1.5, -.9, 1.3, 2.2]])]
    ldr_gates = [0 if lg < 0 else lg for lg in ldr_gates]

    # Plotting
    # --------
    if ax is None:
        ax = plt.figure()
    ax.scatter(log_dna, ldrtxt, c=z, s=10)
    ax.plot([dna_gates[i] for i in [0, 0, 3, 3, 0]],
            [ldr_gates[i] for i in [0, 1, 1, 0, 0]], '-r')
    ax.plot([dna_gates[i] for i in [1, 1, 2, 2, 1]],
            [ldr_gates[i] for i in [0, 1, 1, 0, 0]],
            '-r', linewidth=2)

    ax.plot(g1_g2_pos, [0, 0], 'xk', )
    ax.plot(g1_g2_pos, [0, 0], 'ok', markersize=14, markerfacecolor='None')
    ax.set_xlabel('log10 (DNA content)')
    ax.set_ylabel('LDRtxt intensity')
    dna_lims = get_dnalims(log_dna, x_dna)
    dna_lims = [np.min((dna_lims[0], dna_gates[0]-0.1)),
                np.max((dna_lims[1], dna_gates[3]+0.1))]
    ldr_lims = get_ldrlims(ldrtxt, x_ldr)
    ax.set_xlim(dna_lims)
    ax.set_ylim(ldr_lims)


def live_dead(ldrtxt, ldr_gates,
              dna=None, dna_gates=None,
              x_ldr=None, x_dna=None, ax=None):
    """Assign classification to individual cells as live/dead based on
    ldrtxt and DNA content.
    If ax is not None, plots pie chart of fraction live/dead
    1. alive = selected+others, where selected is within
                 inner DNA gate and within LDR
    2. dead = anything outside of DNA outer gating and LDR gating
    3. total = alive + dead; selected + others + dead

    Parameters
    ----------
    ldrtxt : 1d array
       ldr txt feature across all cells in a well
    ldr_gates : list of floats
    dna : 1d array
       DNA content of cells in a given well
    dna_gates : list of floats
    x_ldr : 1d array
       uniformly distributed 1d grid based on expected
        range of ldr txt
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
    if x_ldr is None:
        mx = np.max(ldrtxt.tolist())+0.01
        x_ldr = np.arange(-0.01, mx, 0.0002)
    outcome = [0] * len(ldrtxt)
    ldr_gates = get_ldrgates(ldrtxt, x_ldr)
    ldr_outer = (ldrtxt < ldr_gates[0]) | (ldrtxt > ldr_gates[1])
    outcome = [-1 if b else 0 for b in ldr_outer]
    dead = np.sum([1 for ot in outcome if ot == -1])
    alive = np.sum([1 for ot in outcome if ot >= 0])
    selected = 'DNA information unavailable'
    others = 'DNA information unavailable'

    if dna is not None:
        log_dna = compute_log_dna(dna, x_dna)
        dna_outermost = (log_dna < dna_gates[0]) | (log_dna > dna_gates[3])
        dna_inner = ((log_dna > dna_gates[1]) &
                     (log_dna < dna_gates[2]) &
                     (ldr_outer == False))
        outcome = [-1 if d else 1 if s else 0
                   for d, s in zip((ldr_outer | dna_outermost), dna_inner)]
        alive = np.sum([1 for ot in outcome if ot >= 0])
        dead = np.sum([1 for s in outcome if s == -1])
        selected = np.sum([1 for s in outcome if s == 1])
        others = np.sum([1 for s in outcome if s == 0])
        if ax is not None:
            ax.pie([selected, others, dead],
                   labels=['selected', 'others', 'dead'],
                   explode=(0.1, 0.1, 0.1), autopct='%1.1f%%')
            ax.axis('equal')
    else:
        if ax is not None:
            ax.pie([alive, dead], labels=['alive', 'dead'],
                   explode=(0.1, 0.1), autopct='%1.1f%%')
            ax.axis('equal')
    return alive, dead, outcome


def plot_summary(ldr, dna, x_ldr=None, well=None):
    """Master summary plot which incorporates above plots as subplots
    
    Parameters
    ----------
    ldr : 1d array
       ldr txt feature across all cells in a well
    dna : 1d array
       DNA content of cells in a given well
    well : str

    Returns
    -------
    fig : matplotlib figure object
    """
    fig = plt.figure()
    gridspec.GridSpec(2, 2)
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=1, rowspan=1)
    ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1, rowspan=1)
    ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=1, rowspan=1)
    ax4 = plt.subplot2grid((2, 2), (1, 1), colspan=1, rowspan=1)
    ldr_gates, ldr_lims = plot_ldr_gating(ldr, x_ldr=x_ldr, ax=ax1)
    dna_gates = get_dna_gating(dna, ldr, ldr_gates, ax=ax2)
    plot_ldr_dna_scatter(dna, ldr, x_ldr=x_ldr, ax=ax3)
    a, d, o = live_dead(ldr, ldr_gates, dna, dna_gates, x_ldr=x_ldr, ax=ax4)
    fig.tight_layout()
    fig.set_size_inches(w=8, h=7)
    if well:
        fig.savefig('dead_cell_filter_%s.png' % well, dpi=300)
    return fig
