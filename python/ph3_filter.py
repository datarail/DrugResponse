import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import smooth
from scipy.stats.mstats import mquantiles as quantile
from findpeaks import get_kde, findpeaks
import matplotlib.gridspec as gridspec


def compute_log_ph3(ph3, x_ph3=None):
    if x_ph3 is None:
        x_ph3 = np.arange(2.5, 8, 0.02)
    ph3_upper_bound = 10 ** x_ph3[-3]
    ph3_lower_bound = 10 ** x_ph3[2]
    ph3_upper_bounded = [d if d < ph3_upper_bound else ph3_upper_bound
                         for d in ph3]
    ph3_bounded = [d if d > ph3_lower_bound else ph3_lower_bound
                   for d in ph3_upper_bounded]
    log_ph3 = np.array([np.log10(d) for d in ph3_bounded])
    return log_ph3



def get_ph3_gates(ph3, cell_identity, x_ph3=None, ph3_cutoff=None):
    if x_ph3 is None:
        x_ph3 = np.arange(2.5, 8, 0.02)
    log_ph3 = compute_log_ph3(ph3, x_ph3)
    if np.any((cell_identity == 1) | (cell_identity == 3)):
        log_ph3_g12 = log_ph3[(cell_identity == 1) | (cell_identity == 3)]
        f_ph3 = get_kde(log_ph3_g12, x_ph3, 4 * (x_ph3[1] - x_ph3[0]))
    else:
        f_ph3 = get_kde(log_ph3, x_ph3, 4 * (x_ph3[1] - x_ph3[0]))
    # if not ph3_cutoff or np.mean(log_ph3 > ph3_cutoff):
    _, peak_loc, peak_width = findpeaks(f_ph3.tolist(), npeaks=3)
    # Enforce that no more than 30% of cells are in M-phase
    min_idx = np.nonzero(np.cumsum(f_ph3)/np.sum(f_ph3) > 0.3)[0][0] - 5
    if np.any(peak_loc > min_idx):
        peak_width = peak_width[np.nonzero(peak_width >= min_idx)[0][0]]
        peak_loc = np.max((peak_loc[np.nonzero(peak_loc >= min_idx)[0][0]],
                           np.nonzero(np.cumsum(f_ph3)/np.sum(f_ph3) > .3)[0][0]
                           ))
    else:
        peak_loc = min_idx
        peak_width = np.max(peak_width)
    
    # find miniminum
    # --------------
    f_ph3_neg = [-x for x in f_ph3[peak_loc:]]
    _, peak_loc_min, _ = findpeaks(f_ph3_neg, npeaks=1)
    peak_loc_min += peak_loc - 1
    ph3_cutoff = x_ph3[math.ceil(np.max((
        np.min((peak_loc_min[0], peak_loc + 9 * peak_width)),
        peak_loc + 2 * peak_width))
                                 )
                    ]
    if not np.any(ph3_cutoff):
        ph3_cutoff = x_ph3[np.nonzero(smooth.smooth(f_ph3, 5, 'flat') >
                                      0.1/len(ph3))[0][0] + 1]
    ph3_lims = (quantile(log_ph3, [5e-3, 0.995]) +
                [(3 * (x_ph3[1] - x_ph3[0])) * q for q in [-1, 10]])
    if np.max(ph3_lims) < ph3_cutoff:
        ph3_lims[1] = ph3_cutoff + 0.02
    return f_ph3, ph3_cutoff, ph3_lims


def evaluate_Mphase(log_ph3, ph3_cutoff, cell_identity, ax=None):
    midx = (log_ph3 >= ph3_cutoff)
    ph3_cell_identity = cell_identity
    ph3_cell_identity[midx] = 4 + ph3_cell_identity[midx]/10
    fractions = {}
    for state, val in zip(['other', 'G1', 'S', 'S_dropout', 'G2', 'M'],
                          [0, 1, 2, 2.1, 3, 4]):
        fractions[state] = np.mean(np.floor(ph3_cell_identity) == (val % 5))
    if ax is None:
        ax = plt.figure()
    ax.pie(fractions.values(), labels=fractions.keys(),  autopct='%1.1f%%')
    ax.axis('equal')
    return fractions


def plot_summary(ph3, cell_identity, x_ph3=None, ph3_cutoff=None, well=None):
    if x_ph3 is None:
        x_ph3 = np.arange(2.5, 8, 0.02)
    log_ph3 = compute_log_ph3(ph3, x_ph3)

    fig = plt.figure()
    gridspec.GridSpec(1, 2)
    ax1 = plt.subplot2grid((1, 2), (0, 0), colspan=1, rowspan=1)
    ax2 = plt.subplot2grid((1, 2), (0, 1), colspan=1, rowspan=1)

    f_ph3, ph3_cutoff, ph3_lims = get_ph3_gates(ph3, cell_identity)
    ax1.plot(x_ph3, np.log10(f_ph3 + np.max(f_ph3)/100) - np.log10(np.max(f_ph3)/100))

    fall = get_kde(log_ph3, x_ph3, 2.5 * (x_ph3[1] - x_ph3[0]))
    ax1.plot(x_ph3, np.log10(fall + np.max(fall)/100) - np.log10(np.max(fall)/100), '--')

    ax1.plot([ph3_cutoff, ph3_cutoff], [0, np.log10(np.max(f_ph3))], '-')
                

    ax1.set_xlim(ph3_lims)
    ax1.set_ylim([0, np.log10(np.max(f_ph3)) - np.log10(np.max(f_ph3)/100)+0.1])
    ax1.set_xlabel('log10 (pH3)')
    ax1.set_ylabel('kernel density estimate')

    fractions = evaluate_Mphase(log_ph3, ph3_cutoff, cell_identity, ax=ax2)
    if well:
        fig.savefig('pH3_%s.png' % well, dpi=300)
    return fractions
