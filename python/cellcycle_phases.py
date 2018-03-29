from findpeaks import findpeaks, get_kde
import numpy as np
from scipy.stats.mstats import mquantiles as quantile
import matplotlib.pyplot as plt
import smooth
from scipy.stats import gaussian_kde
import math
import pandas as pd
import accum
from scipy.ndimage.filters import maximum_filter
from scipy.spatial.distance import pdist, squareform
from scipy.stats import norm
import matplotlib.gridspec as gridspec


df = pd.read_table('training_sample_object_level.txt')
edu = df['Nuclei Selected - EdUINT'].tolist()
edu = np.array(edu)
dna = df['Nuclei Selected - DNAcontent'].tolist()
dna = np.array(dna)


# Gating based on EdU
# -------------------


def get_edu_gates(edu, px_edu, plotting=False, ax=None):
    """ Returns estimate of max EdU for G1 gating and min EdU for S phase gating
    """
    #if not px_edu:
    #    px_edu = np.arange(-0.2, 5.3, .02)
    x_edu = np.arange(-200, 4e3+1, 1)
    # Note: Bandwidth = 90 reproduced MATLAB output
    f_edu = get_kde(edu, x_edu, bandwidth=90)
    peak_amp, peak_loc, peak_width = findpeaks(f_edu.tolist(), npeaks=2)
    peak_amp = peak_amp[peak_amp > np.max(f_edu)/10]
    peak_loc = peak_loc[peak_amp > np.max(f_edu)/10]
    peak_width = peak_width[peak_amp > np.max(f_edu)/10]
    if peak_loc.size == 0:
        x_edu = np.arange(-200, 2e4+1, 1)
        f_edu = get_kde(edu, x_edu, bandwidth=90)
        peak_amp, peak_loc, peak_width = findpeaks(f_edu.tolist(), npeaks=2)
        peak_amp = peak_amp[peak_amp > np.max(f_edu)/10]
        peak_loc = peak_loc[peak_amp > np.max(f_edu)/10]
        peak_width = peak_width[peak_amp > np.max(f_edu)/10]
    peak_width = peak_width[np.argmin(peak_loc)]
    peak_loc = x_edu[math.ceil(np.min(peak_loc))]

    # Find location of minimum on right
    if np.any(edu > (peak_loc + 30)):
        edu_higher = edu[edu > peak_loc + 30]
        f2_edu = get_kde(edu_higher, x_edu, bandwidth=510)
        f2_edu_neg = [-x for x in f2_edu]
        _, peak_trough, _ = findpeaks(f2_edu_neg, npeaks=2)
        peak_trough = x_edu[math.ceil(
            peak_trough[np.argmin(np.abs([x - 500 for x in peak_trough]))])]
        peak_trough = np.max([peak_trough, peak_loc+3*peak_width])
    else:
        peak_trough = peak_loc + 3 * peak_width
    # Edu offset
    # ** Not entirely clear to me yet
    offset_edu = np.max((peak_loc-1.5 * peak_width, 1))

    log_edu = compute_log_edu(edu, px_edu, offset_edu)

    # EdU max for G1 gating
    edu_g1_max = np.max((
        np.log10(peak_trough - offset_edu),  # Expected EdU max for G1 (optn 1)
        quantile(log_edu, 0.2) + 0.1  # Expected EdU max  for G1 (optn 2)
    ))

    # Edu  min for S phase gating
    edu_s_min = np.max((
        np.log10(peak_loc + 2 * peak_width - offset_edu),  # Expected EdU min for S phase (optn 1)
        edu_g1_max - 0.1  # Expected Edu min for S phase (option 2)
        ))

    # Expected differene between G1/G2 and S
    edu_shift = np.max((
        (np.log10(peak_loc + 2 * peak_width - offset_edu) -  # Expected EdU min for S phase (optn 1)
         np.log10(np.max((peak_loc - offset_edu, 1)))),  # G1 peak loc (offset)
        1
        ))

    # Plotting
    # --------
    if plotting:
        if ax is None:
            ax = plt.figure()            
        idx = np.random.permutation(len(edu))
        idx = idx[:np.min((len(edu), 1000))]
        edu = np.array(edu)
        ax.plot(edu[idx], log_edu[idx], '.c')
        px_edu = np.arange(-0.2, 5.3, .02)
        ax.plot(x_edu, np.max(px_edu) * (f_edu/np.max(f_edu)), 'k-')
        ax.plot(x_edu, np.max(px_edu) * (f2_edu/np.max(f2_edu)), 'k--')
        ax.plot([-200, 300, np.nan, -100, 500],
                 [edu_s_min, edu_s_min, np.nan,
                  edu_shift+np.log10(np.max((peak_loc-offset_edu, 1))),
                  edu_shift+np.log10(np.max((peak_loc-offset_edu, 1)))], '-r')
        ax.plot([offset_edu, offset_edu], [0, 5], ':r')
        ax.plot([100, peak_trough, peak_trough],
                 [edu_g1_max, edu_g1_max, 0], '--r')
        ax.set_ylim((px_edu[0], px_edu[-1]))
        ax.set_xlim([-200, np.max((peak_loc + 5 * peak_width, 500))])
        ax.set_xlabel('EdU intensity')
        ax.set_ylabel('log10 (EdU')
    return edu_shift, offset_edu, edu_g1_max, edu_s_min


def compute_log_dna(dna, x_dna=None):
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


def compute_log_edu(edu, px_edu, offset_edu):
    edu_upper_bound = 10 ** px_edu[-3]
    edu_lower_bound = 10 ** px_edu[2]
    edu_offsetted = [e - offset_edu for e in edu]
    edu_upper_bounded = [e if e < edu_upper_bound else edu_upper_bound
                         for e in edu_offsetted]
    edu_bounded = [e if e > edu_lower_bound else edu_lower_bound
                   for e in edu_upper_bounded]
    log_edu = np.array([np.log10(e) for e in edu_bounded])
    return log_edu


def plot_edu_dna_scatter(dna, edu, offset_edu,
                         dna_gates, edu_gates,
                         dna_lims, edu_lims,
                         x_dna=None, px_edu=None,
                         ax=None):
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    if px_edu is None:
        px_edu = np.arange(-0.2, 5.3, .02)
    log_dna = compute_log_dna(dna, x_dna)
    log_edu = compute_log_edu(edu, px_edu, offset_edu)
    xy = np.vstack([log_dna, log_edu])
    z = gaussian_kde(xy)(xy)
    #scatter_fig = plt.figure()
    ax.scatter(log_dna, log_edu, c=z, s=10)
    ax.plot([dna_gates[i] for i in [0, 0, 3, 3, 0, 0, 1, 1, 0, 2, 2, 3]],
             [-1, edu_gates[1], edu_gates[1], -1, np.nan, edu_gates[0], edu_gates[0],
              -1, np.nan, -1, edu_gates[0], edu_gates[0]], '--', color='red')
    ax.set_xlabel('log10 (DNA content)')
    ax.set_ylabel('log10 (EdU)')
    ax.set_xlim(dna_lims)
    ax.set_ylim(edu_lims)
    # return scatter_fig

# Find peaks in 2-dimensions (EdU and DNA)
# ----------------------------------------


def histc(X, bins):
    """ Replicate MATLAB histc function that returns counts per bin and
    bin index for each data point
    """
    map_to_bins = np.digitize(X, bins)
    r = np.zeros(bins.shape)
    for i in map_to_bins:
        r[i-1] += 1
    return [r, map_to_bins]


def smooth_1d(y, lm=5):
    m, n = y.shape
    e = np.identity(m)
    d1 = np.diff(e, 1).T
    d2 = np.diff(d1.T, 1).T
    p = (lm ** 2) * np.matmul(d2.T, d2) + 2 * lm * np.matmul(d1.T, d1)
    z = np.linalg.lstsq((e + p), y)
    return z[0]


def imregionalmax(f):
    conn_8 = np.ones((8, 8))
    regional_max = maximum_filter(f, footprint=conn_8)
    peak_2d = (f == regional_max)
    return peak_2d


def get_2d_histogram(log_dna, x_dna, log_edu, px_edu):
    """ Count the log intensity in each 2-D bin
    """
    # ---------------------------------------
    # get bin index for each log intensity value
    _, bin_indeces_dna = histc(log_dna, x_dna)
    _, bin_indeces_edu = histc(log_edu, px_edu)

    bin_indeces = np.array([bin_indeces_edu, bin_indeces_dna]).T
    vals = np.array([1] * bin_indeces.shape[0])
    h = accum.accum(bin_indeces, vals, size=[len(px_edu), len(x_dna)])
    h[h == 1] = 0
    h = h/len(log_dna)
    return h


def get_2D_peak(h, x_dna, px_edu, nsmooth=5):
    """ Return peaks candidates from 2D readout of DNA and EdU
    """
    g = smooth_1d(h, nsmooth)
    f = smooth_1d(g.T, nsmooth).T
    peak_2d = imregionalmax(f)
    x, y = np.nonzero(peak_2d)
    
    pre_peak_candidates = np.array([x_dna[y], px_edu[x], f[peak_2d]]).T
    peak_candidates = pre_peak_candidates[
        ((px_edu[x] > (nsmooth + 2) * (px_edu[1] - px_edu[0])) &
         (f[peak_2d] > 1e-5) & (f[peak_2d] > np.max(f[peak_2d]/30))),
        :]

    # Sort by 3rd column (descending order)
    peak_candidates = peak_candidates[peak_candidates[:, 2].argsort()[::-1]]
    return peak_candidates, pre_peak_candidates.shape[0]


def iterate_2D_peak(h, x_dna, px_edu, nsmooth=5):
    peak_candidates, lenp = get_2D_peak(h, x_dna, px_edu, nsmooth)
    if peak_candidates.shape[0] < 2:
        nsmooth = 0.5 * nsmooth
        peak_candidates, lenp = get_2D_peak(h, x_dna, px_edu, nsmooth)
    elif (lenp > 2 * peak_candidates.shape[0]) | (lenp > 5):
        nsmooth = 2 * nsmooth
        peak_candidates, lenp = get_2D_peak(h, x_dna, px_edu, nsmooth)
    return peak_candidates


def plot_2D_peaks(log_dna, x_dna, edu, px_edu,
                  f, peak_candidates, phase_candidates,
                  dna_gates, dna_lims,
                  edu_gates, edu_lims,
                  nsmooth=5, ax=None):
    # edu_shift, offset_edu, edu_g1_max, edu_s_min = get_edu_gates(edu, px_edu)
    # log_edu = compute_log_edu(edu, px_edu, offset_edu)
    # h = get_2d_histogram(log_dna, x_dna, log_edu, px_edu)
    # g = smooth_1d(h, nsmooth)
    # f = smooth_1d(g.T, nsmooth).T
    # peak_candidates = iterate_2D_peak(h, x_dna, px_edu, nsmooth)
    # edu_shift, offset_edu, edu_g1_max, edu_s_min = get_edu_gates(edu, px_edu)
    # phase_candidates = get_phase_candidates(peak_candidates,
    #                                        edu_shift, edu_s_min)
    if ax is None:
        ax = plt.figure()
    ax.pcolor(x_dna, px_edu, f)
    ax.plot(peak_candidates[:, 0], peak_candidates[:, 1], 'ok',
            markersize=4, markerfacecolor='None')
    phases = ['G1', 'S', 'G2']
    for i, phase in enumerate(phases):
        ax.text(phase_candidates[i, 0], phase_candidates[i, 1], phase,
                fontsize=14, color='red')
    ax.plot([dna_gates[i] for i in [0, 0, 3, 3, 0, 0, 1, 1, 0, 2, 2, 3]],
            [-1, edu_gates[1], edu_gates[1], -1, np.nan,
             edu_gates[0], edu_gates[0],
             -1, np.nan, -1, edu_gates[0], edu_gates[0]],
            linewidth=2, color='red')
    ax.set_xlim(dna_lims)
    ax.set_ylim(edu_lims)
    ax.set_xlabel('log (DNA content)')
    ax.set_ylabel('log (EdU intensity)')
    # ax.set_xlim(quantile(x_dna, [0.25, 0.75]))
    # ax.set_ylim(quantile(px_edu, [0.1, 0.75]))


# Assign cell cycle phase based on peak candidates
# ------------------------------------------------
def get_phase_candidates(peak_candidates, edu_shift, edu_s_min):

    phase_candidates = np.zeros((3, 2))
    edu_peak_bool = (((peak_candidates[:, 1] - np.min(peak_candidates[:, 1]))
                      > edu_shift) &
                     (peak_candidates[:, 1] > (edu_s_min + 0.2 * edu_shift)))
    if np.any(edu_peak_bool):
        repmat1 = np.matlib.repmat(
            peak_candidates[:, 0].reshape(peak_candidates.shape[0], 1),
            1, peak_candidates.shape[0])
        repmat2 = np.matlib.repmat(
            peak_candidates[:, 2].reshape(peak_candidates.shape[0], 1),
            1, peak_candidates.shape[0])
        s_phase_candidates = peak_candidates[np.argmax(
            edu_peak_bool *
            np.sum(np.abs(repmat1 - repmat1.T < np.log10(2) * 0.75) *
                   repmat2, axis=0).T +
            peak_candidates[:, 2]
        ), [0, 1]
        ]
        phase_candidates[1, :] = s_phase_candidates
        g1_candidates = peak_candidates[((peak_candidates[:, 0] < s_phase_candidates[0] + np.log10(2) * 0.05) &
                                        (peak_candidates[:, 0] > s_phase_candidates[0] - 0.75*np.log10(2)) &
                                         (peak_candidates[:, 1] < s_phase_candidates[1] - edu_shift)), :]
        if np.any(g1_candidates):
            phase_candidates[0, :] = g1_candidates[np.argmax(g1_candidates[:, 2]), [0, 1]]

        g2_candidates = peak_candidates[((peak_candidates[:, 0] > np.nanmean(phase_candidates[:2, 0])) &
                                         (peak_candidates[:, 0] < np.nanmean(phase_candidates[:2, 0]) + np.log10(2)) &
                                         (peak_candidates[:, 1] < phase_candidates[1, 1] - edu_shift)), :]
        if np.any(g2_candidates):
            phase_candidates[2, :] = g2_candidates[np.argmax(g2_candidates[:, 2]), [0, 1]]
    else:
        # Most likely no S-phase_candidates, therefore assign 2 highest peaks as G1 ans G2 based on DNA
        if np.any(peak_candidates):
            if peak_candidates.shape[0] == 1:  # If only 1 peak, just assign to G1
                phase_candidates[0, :] = peak_candidates[0, :2]
            else: # take the ones that are best seperated
                p1 = peak_candidates[:, 0]
                p2 = np.zeros(peak_candidates.shape[0])
                p3 = np.concatenate((p1, p2)).reshape(2, len(p2)).T
                # pdist(p3)
                repmat3 = np.matlib.repmat(
                    peak_candidates[:, 2].reshape(peak_candidates.shape[0], 1),
                    1, peak_candidates.shape[0])
                repmat4 = np.matlib.repmat(
                    peak_candidates[:, 2].reshape(peak_candidates.shape[0], 1).T,
                    peak_candidates.shape[0], 1)
                pk_dist = (squareform(pdist(p3)) > 0.6*np.log10(2)) * (repmat3 + repmat4)
                pk_dist_idx1, pk_dist_idx2 = (np.nonzero(pk_dist == np.max(pk_dist)))[0]
                if peak_candidates[pk_dist_idx1, 0] > peak_candidates[pk_dist_idx2, 0]:
                    phase_candidates[0, :] = peak_candidates[pk_dist_idx2, [0, 1]]
                    phase_candidates[2, :] = peak_candidates[pk_dist_idx1, [0, 1]]
                else:
                    phase_candidates[0, :] = peak_candidates[pk_dist_idx1, [0, 1]]
                    phase_candidates[2, :] = peak_candidates[pk_dist_idx2, [0, 1]]
    return phase_candidates


# Working with each channel sequentially
# --------------------------------------
def get_g1_dna_peak(log_dna, x_dna, log_edu, edu_shift,
                    edu_s_min, edu_g1_max, phase_candidates,
                    ax=None):
    f_dna = get_kde(log_dna, x_dna)
    log_dna_low_edu = log_dna[(log_edu < edu_s_min+0.2*edu_shift) &
                              (log_edu < edu_g1_max)]
    f_dna_low_edu = get_kde(log_dna_low_edu, x_dna)
    peak_amp, peak_loc, _ = findpeaks(f_dna_low_edu.tolist())
    peak_loc = peak_loc[peak_amp > np.max(peak_amp/10)]
    dna_g1_loc = x_dna[peak_loc[:3]]

    if len(dna_g1_loc) > 1:
        if phase_candidates[0, 0]:
            dna_g1_loc = dna_g1_loc[np.argmin(np.abs(
                dna_g1_loc - phase_candidates[0, 0]))]
        elif phase_candidates[1, 0]:
            dna_g1_loc = np.max(dna_g1_loc[
                dna_g1_loc > phase_candidates[1, 0]])
        else:
            dna_g1_loc = np.min(dna_g1_loc)
    if not np.any(dna_g1_loc):
        dna_g1_loc = np.nanmin(phase_candidates[:, 0] - np.log10(1.2))
    if ax is None:    
        ax = plt.figure('dna_figure')
    ax.plot(x_dna, f_dna)
    ax.plot(dna_g1_loc, .1, 'xk')
    ax.set_xlabel('log (DNA content)')
    ax.set_ylabel('kernel density estimate')
    return dna_g1_loc# , dna_fig


# Working with EdU channel
# ------------------------


def get_low_edu_peaks(log_edu, px_edu, edu_shift,
                      edu_g1_max,
                      log_dna, dna_g1_loc, nsmooth=5):
    f_edu = get_kde(log_edu, px_edu)
    log_edu_low_bool = ((log_dna > dna_g1_loc - 1) &
                        (log_dna < dna_g1_loc + 0.1) &
                        (log_edu > 2 * nsmooth * (px_edu[1] - px_edu[0])) &
                        (log_edu < edu_g1_max))
    if not np.any(log_edu_low_bool):
        log_edu_low_bool = ((log_dna > dna_g1_loc - 1) &
                            (log_dna < dna_g1_loc + 0.1) &
                            (log_edu < edu_g1_max))
    f_edu_low = get_kde(log_edu[log_edu_low_bool], px_edu)
    bin_counts, _ = histc(log_edu[log_edu_low_bool], px_edu)
    # Check discrepency in array length when using 3
    f_edu_low[[True, smooth.smooth(bin_counts, 2.99, 'flat') <= 1/3]] = 0

    edu_amp, edu_loc, _ = findpeaks(smooth.smooth(f_edu_low, nsmooth).tolist(),
                                    npeaks=2)
    if np.any(edu_loc):
        edu_loc = edu_loc[edu_amp > 0.3 * np.max(edu_amp)]
        low_edu_peaks = px_edu[edu_loc[np.argmin(edu_loc)]]
    else:
        low_edu_peaks = np.median(log_edu[log_edu_low_bool])
    return low_edu_peaks


# Get peaks with high EdU values
# ------------------------------
def get_high_edu_peaks(log_edu, px_edu, edu_shift,
                       low_edu_peaks, log_dna, dna_g1_loc,
                       nsmooth=5):
    high_edu_bool = ((log_dna > dna_g1_loc - np.log10(2)/2) &
                     (log_dna < dna_g1_loc + np.log10(2) * 1.5) &
                     (log_edu > low_edu_peaks + edu_shift * 0.8))
    if np.any(high_edu_bool):
        f_edu_high = get_kde(log_edu[high_edu_bool], px_edu)
        edu_amp, edu_loc, _ = findpeaks(smooth.smooth(
            f_edu_high, nsmooth, 'flat').tolist())
        # Remove lesser peaks
        high_edu_loc = edu_loc[edu_amp > np.max(edu_amp/10)]

        high_edu_bool = px_edu[high_edu_loc] > low_edu_peaks + edu_shift
        if np.any(high_edu_bool):
            high_edu_peaks = px_edu[high_edu_loc[
                np.nonzero(high_edu_bool)[0][0]]]
        else:
            high_edu_peaks = low_edu_peaks + edu_shift

        f_edu = get_kde(log_edu, px_edu)
        neg_f_edu = [-x for x in f_edu]
        _, edu_loc, _ = findpeaks(neg_f_edu)
        edu_cutoff = px_edu[edu_loc[
            ((np.nonzero(px_edu[edu_loc] > low_edu_peaks) &
              (px_edu[edu_loc] < high_edu_peaks)))[0][0]
        ]]
        if not np.any(edu_cutoff):
            edu_cutoff = px_edu[np.argmin(smooth.smooth(f_edu.T, nsmooth, 'flat').T +
                                          ((px_edu < low_edu_peaks) |
                                           (px_edu > high_edu_peaks)))]
        
        edu_lims = [px_edu[2],
                    np.min((2 * high_edu_peaks - edu_cutoff, px_edu[-2]))]
        # dna_s_loc = get_s_phase_dna_peaks(log_dna, x_dna, dna_g1_loc,
        #                                   log_edu, edu_cutoff)
    else:
        high_edu_peaks = low_edu_peaks + edu_shift
        edu_cutoff = np.mean(low_edu_peaks, high_edu_peaks)
        edu_lims = [-0.02,
                    np.min((2*high_edu_peaks - edu_cutoff + 0.1, px_edu[-2]))]
        # dna_s_loc = dna_g1_loc + 0.5 * np.log10(2)
    edu_gates = [edu_cutoff,
                 high_edu_peaks + np.max((high_edu_peaks-edu_cutoff, 1))]
    edu_lims[1] = np.max((edu_lims[1], edu_gates[1]+0.1))
    edu_peaks = [low_edu_peaks, high_edu_peaks]
    return edu_peaks, edu_cutoff, edu_lims, edu_gates


def get_s_phase_dna_loc(log_dna, x_dna,
                        dna_g1_loc,  # dna_fig,
                        log_edu, edu_cutoff,
                        nsmooth=5, ax=None):
    high_dna_bool = ((log_dna > dna_g1_loc - np.log10(2) * 0.5) &
                     (log_dna < dna_g1_loc + np.log10(2) * 1.5) &
                     (log_edu > edu_cutoff))
    if np.any(high_dna_bool):
        f_dna = get_kde(log_dna[high_dna_bool], x_dna, bandwidth=0.0317)
        dna_amp, dna_loc, _ = findpeaks(
            smooth.smooth(f_dna, nsmooth, 'flat').tolist(),
            npeaks=1)
        dna_s_loc = x_dna[dna_loc[0]]
    else:
        dna_s_loc = dna_g1_loc + 0.5 * np.log10(2)
    # ax.figure(dna_fig.number)
    ax.plot(x_dna, f_dna, '-.')
    return dna_s_loc #, dna_fig


def get_g2_dna_loc(log_dna, x_dna, log_edu, edu_cutoff,
                   dna_g1_loc,  # dna_fig,
                   phase_candidates, nsmooth, ax):
    high_dna_bool = ((log_dna > dna_g1_loc + 0.4 * np.log10(2)) &
                     (log_edu < edu_cutoff))
    f_dna = get_kde(log_dna[high_dna_bool], x_dna)

    peak_amp, peak_loc, _ = findpeaks(smooth.smooth(f_dna, nsmooth).tolist())
    peak_loc = peak_loc[peak_amp > np.max(peak_loc/10)]
    g2_loc_candidates = x_dna[peak_loc]
    g2_loc_candidates = g2_loc_candidates[g2_loc_candidates >
                                          (dna_g1_loc + 0.5*np.log10(2))]
    if len(g2_loc_candidates) > 1:
        if np.any(phase_candidates[2, 0]):
            g2_loc = g2_loc_candidates[np.argmin(
                np.abs(g2_loc_candidates - phase_candidates[2, 0]))]
        else:
            bc = ((np.any(phase_candidates[1, 0])) &
                  (np.any(g2_loc_candidates > phase_candidates[1, 0])))
            if bc:
                g2_loc = g2_loc_candidates[
                    g2_loc_candidates > phase_candidates[1, 0]]
            g2_loc = g2_loc[np.argmin(np.abs(
                g2_loc - dna_g1_loc - np.log10(2)))]
    elif len(g2_loc_candidates) == 1:
                g2_loc = g2_loc_candidates
    else:
            g2_loc = dna_g1_loc + np.log10(2)
    # plt.figure(dna_fig.number)
    ax.plot(x_dna, f_dna, ':')
    return g2_loc#, dna_fig


def get_dna_cutoff(log_dna, x_dna, log_edu, edu_cutoff,
                   dna_g1_loc, dna_s_loc,# dna_fig,
                   phase_candidates, nsmooth, ax):
    high_dna_bool = ((log_dna > dna_g1_loc + 0.4 * np.log10(2)) &
                     (log_edu < edu_cutoff))
    if np.any(high_dna_bool):
        dna_g2_loc = get_g2_dna_loc(log_dna, x_dna, log_edu, edu_cutoff,
                                    dna_g1_loc,  # dna_fig,
                                    phase_candidates, nsmooth, ax=ax)
        f_dna = get_kde(log_dna[log_edu < edu_cutoff], x_dna)
        smooth_f_dna = smooth.smooth(f_dna, nsmooth, 'flat')
        _, peak_loc, _ = findpeaks([-x for x in smooth_f_dna])
        dna_cutoff = x_dna[peak_loc[((x_dna[peak_loc] > dna_g1_loc) &
                                     (x_dna[peak_loc] < dna_g2_loc))]]
        
        if not np.any(dna_cutoff):
            dna_cutoff = np.min((np.max((dna_s_loc, dna_g1_loc + 0.02)),
                                 dna_g2_loc - 0.02))
        elif len(dna_cutoff) >= 1:
            dna_cutoff = dna_cutoff[0]
    else:
        dna_cutoff = dna_g1_loc + 0.3 * np.log10(2)
        dna_g2_loc = dna_g1_loc + np.log10(2)
    # plt.figure(dna_fig.number)
    ax.plot([dna_g1_loc, dna_s_loc, dna_g2_loc], [0, 0, 0], 'xk')
    ax.plot(dna_cutoff, 0, 'xk')
    return dna_cutoff, dna_g2_loc#, dna_fig


def get_normal_dist(data):
    mu, std = norm.fit(data)
    xmin, xmax = data.min(), data.max()
    x = np.linspace(xmin, xmax, len(data))
    p = norm.pdf(x, mu, std)
    return p, mu, std


# Find cells dropping in S-phase
# ------------------------------
def get_dna_gates(log_dna, x_dna, dna_g1_loc, dna_g2_loc, dna_cutoff,
                  log_edu, edu_cutoff):
    hg1 = np.abs(((log_dna - dna_g1_loc < 0.3 * np.log10(2)) &
                  (log_edu < edu_cutoff)))
    if np.sum(hg1) > 10:
        norm_fit_g1, mu, std = get_normal_dist(log_dna[hg1])
        g1_min_width = np.max((norm.ppf(0.9, mu, std) - norm.ppf(0.1, mu, std),
                               0.05))
        g1_lim = np.min((norm.ppf(0.99, mu, std),
                         dna_cutoff - 0.1 * np.log10(2)))
    else:
        g1_min_width = 0.05
        g1_lim = np.min((dna_cutoff - 0.1 * np.log10(2),
                         np.mean((dna_cutoff, dna_g1_loc))))

    hg2 = np.abs(((log_dna - dna_g2_loc < 0.3 * np.log10(2)) &
                  (log_edu < edu_cutoff)))
    if np.sum(hg2) > 10:
        norm_fit_g2, mu, std = get_normal_dist(log_dna[hg2])
        g2_min_width = np.max((norm.ppf(0.9, mu, std) - norm.ppf(0.1, mu, std),
                               0.05))
        g2_lim = np.max((norm.ppf(0.1, mu, std),
                         dna_cutoff + 0.1 * np.log10(2)))
    else:
        g2_min_width = 0.05
        g2_lim = np.min((dna_cutoff + 0.1 * np.log10(2),
                         np.mean((dna_cutoff, dna_g2_loc))))

    d1 = dna_cutoff - dna_g1_loc
    d2 = dna_g2_loc - dna_cutoff
    dna_lims = np.array([np.max((dna_g1_loc - 3 * d1, x_dna[1])),
                         np.min((dna_g2_loc + 3 * d2, x_dna[-2]))])
    

    dna_gates = np.array([dna_g1_loc - d1, g1_lim, g2_lim, dna_g2_loc+d2])
    
    if dna_gates[1] - dna_gates[0] < g1_min_width:
        dna_gates[:2] = [x * g1_min_width + np.mean(dna_gates[:2])
                         for x in [-0.6, -.4]]
    if dna_gates[3] - dna_gates[2] < g1_min_width:  # CHECK ########
        dna_gates[2:] = [x * g2_min_width + np.mean(dna_gates[2:])
                         for x in [-0.4, .6]]
    if dna_gates[2] < dna_gates[1]:
        dna_gates[1:3] = np.mean(dna_gates[1:3])

    dna_gl = list(dna_lims) + [g-0.1 for g in dna_gates]
    dna_gl2 = list(dna_lims) + [g+0.1 for g in dna_gates]
    dna_lims = [np.min(dna_gl), np.max(dna_gl2)]
    
    return dna_gates, dna_lims


def evaluate_cell_cycle_phase(log_dna, dna_gates, x_dna, dna_peaks,
                              log_edu, edu_gates, px_edu, edu_peaks,
                              nsmooth=5, ax=None):
    cell_id = (1 * ((log_dna > dna_gates[0]) &  # G1
                    (log_dna < dna_gates[1]) &
                    (log_edu < edu_gates[0])) +
               2 * ((log_dna >= dna_gates[0]) &  # S
                    (log_dna < dna_gates[3]) &
                    (log_edu >= edu_gates[0]) &
                    (log_edu < edu_gates[1])) +
               2.1 * ((log_dna >= dna_gates[1])  # S dropout
                      & (log_dna < dna_gates[2]) &
                      (log_edu < edu_gates[0])) +
               3 * ((log_dna >= dna_gates[2]) &  # G2
                    (log_dna < dna_gates[3]) &
                    (log_edu < edu_gates[0])))
    fractions = {}
    for state, val in zip(['other', 'G1', 'S', 'S_dropout', 'G2'],
                          [0, 1, 2, 2.1, 3]):
        fractions[state] = np.mean(cell_id == (val % 4))
    
    for ig in np.arange(1, 4):
        if sum(cell_id == ig) > 10:
            
            f_dna = get_kde(log_dna[cell_id == ig], x_dna)
            _, dna_loc, _ = findpeaks(
                smooth.smooth(f_dna, 3 * nsmooth).tolist(),
                npeaks=1)
            dna_peaks[ig-1] = x_dna[dna_loc]
            
            f_edu = get_kde(log_edu[cell_id == ig], px_edu)
            _, edu_loc, _ = findpeaks(
                smooth.smooth(f_edu, 3 * nsmooth).tolist(),
                npeaks=1)
            edu_peaks[ig-1] = px_edu[edu_loc]
            
        else:
            
            dna_peaks[ig-1] = np.mean(dna_gates[ig-1:2])
            edu_peaks[ig-1] = np.mean((edu_gates[0], (ig == 2)*edu_gates[1]))
        edu_peaks[1] = np.max((edu_peaks[1], edu_gates[0] + 0.1))
    peaks = [dna_peaks, edu_peaks]
    if ax is None:
        ax = plt.figure()
    ax.pie(fractions.values(), labels=fractions.keys(),  autopct='%1.1f%%')
    ax.axis('equal')
    return fractions, cell_id, peaks


def plot_summary(dna, edu, x_dna=None, px_edu=None, well=None):
    fig = plt.figure()
    gridspec.GridSpec(2, 3)
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=1, rowspan=1)
    ax2 = plt.subplot2grid((2, 3), (0, 1), colspan=1, rowspan=1)
    ax3 = plt.subplot2grid((2, 3), (0, 2), colspan=1, rowspan=1)
    ax4 = plt.subplot2grid((2, 3), (1, 0), colspan=2, rowspan=1)
    ax5 = plt.subplot2grid((2, 3), (1, 2), colspan=1, rowspan=1)
    
    if px_edu is None:
        px_edu = np.arange(-0.2, 5.3, .02)
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    # if ax is None:
    #     ax = plt.figure()
    log_dna = compute_log_dna(dna, x_dna)

    edu_shift, offset_edu, edu_g1_max, edu_s_min = get_edu_gates(edu, px_edu,
                                                                 plotting=True, ax=ax1)
    log_edu = compute_log_edu(edu, px_edu, offset_edu)
    h = get_2d_histogram(log_dna, x_dna, log_edu, px_edu)
    g = smooth_1d(h, 5)
    f = smooth_1d(g.T, 5).T

    peak_candidates = iterate_2D_peak(h, x_dna, px_edu, nsmooth=5)
    phase_candidates = get_phase_candidates(peak_candidates,
                                            edu_shift, edu_s_min)

    dna_g1_loc = get_g1_dna_peak(log_dna, x_dna, log_edu, edu_shift,
                                 edu_s_min, edu_g1_max,
                                 phase_candidates, ax=ax2)

    low_edu_peaks = get_low_edu_peaks(log_edu, px_edu, edu_shift,
                                      edu_g1_max,
                                      log_dna, dna_g1_loc, nsmooth=5)
    edu_peaks, edu_cutoff, edu_lims, edu_gates = get_high_edu_peaks(log_edu, px_edu,
                                                                    edu_shift,
                                                                    low_edu_peaks,
                                                                    log_dna, dna_g1_loc,
                                                                    nsmooth=5)

    dna_s_loc = get_s_phase_dna_loc(log_dna, x_dna, dna_g1_loc,
                                    log_edu, edu_cutoff, ax=ax2)

    dna_cutoff, dna_g2_loc = get_dna_cutoff(log_dna, x_dna, log_edu,
                                            edu_cutoff, dna_g1_loc,
                                            dna_s_loc,  # dna_fig,
                                            phase_candidates, 5, ax=ax2)
    dna_gates, dna_lims = get_dna_gates(log_dna, x_dna, dna_g1_loc, dna_g2_loc,
                                        dna_cutoff, log_edu, edu_cutoff)
    ax2.plot([dna_gates[i] for i in [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3]],
             [np.max(ax2.get_ylim()) * i for i in [0, 1, np.nan, 0, 1, np.nan,
                                               0, 1, np.nan, 0, 1]],
             '--', color='red')
    ax2.set_xlim(dna_lims)

    dna_peaks = [dna_g1_loc, dna_s_loc, dna_g2_loc]
    edu_peaks = edu_peaks + [edu_peaks[0]]

    plot_2D_peaks(log_dna, x_dna, edu, px_edu,
                  f, peak_candidates, phase_candidates,
                  dna_gates, dna_lims,
                  edu_gates, edu_lims,
                  nsmooth=5, ax=ax3)
    
    plot_edu_dna_scatter(dna, edu, offset_edu,
                         dna_gates, edu_gates,
                         dna_lims, edu_lims,
                         x_dna=None, px_edu=None,
                         ax=ax4)

    fractions, cell_id, peaks = evaluate_cell_cycle_phase(log_dna, dna_gates,
                                                          x_dna, dna_peaks,
                                                          log_edu, edu_gates,
                                                          px_edu, edu_peaks,
                                                          nsmooth=5, ax=ax5)
    fig.tight_layout()
    fig.set_size_inches(w=10, h=6)
    if well:
        fig.savefig('cell_cycle_phases_%s.png' % well, dpi=300)
    return fractions, cell_id
