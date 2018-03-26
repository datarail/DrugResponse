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


df = pd.read_table('training_sample_object_level.txt')
edu = df['Nuclei Selected - EdUINT'].tolist()
edu = np.array(edu)
dna = df['Nuclei Selected - DNAcontent'].tolist()
dna = np.array(dna)


# Gating based on EdU
# -------------------


def get_edu_gates(edu, px_edu=None):
    """ Returns estimate of max EdU for G1 gating and min EdU for S phase gating
    """
    if not px_edu:
        px_edu = np.arange(-0.2, 5.3, .02)
    x_edu = np.arange(-200, 4e3, 1)
    # Note: Bandwidth = 100 reproduced MATLAB output
    f_edu = get_kde(edu, x_edu, bandwidth=100)
    peak_amp, peak_loc, peak_width = findpeaks(f_edu.tolist(), npeaks=2)
    peak_amp = peak_amp[peak_amp > np.max(f_edu)/10]
    peak_loc = peak_loc[peak_amp > np.max(f_edu)/10]
    peak_width = peak_width[peak_amp > np.max(f_edu)/10]
    if peak_loc.size == 0:
        x_edu = np.arange(-200, 2e4, 1)
        f_edu = get_kde(edu, x_edu, bandwidth=100)
        peak_amp, peak_loc, peak_width = findpeaks(f_edu.tolist(), npeaks=2)
        peak_amp = peak_amp[peak_amp > np.max(f_edu)/10]
        peak_loc = peak_loc[peak_amp > np.max(f_edu)/10]
        peak_width = peak_width[peak_amp > np.max(f_edu)/10]
    peak_width = peak_width[np.argmin(peak_loc)]
    peak_loc = x_edu[math.ceil(np.min(peak_loc))]

    # Find location of minimum on right
    if np.any(edu > (peak_loc + 30)):
        edu_higher = np.array([e for e in edu if e > (peak_loc + 30)])
        f2_edu = get_kde(edu_higher, x_edu)
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
    idx = np.random.permutation(len(edu))
    idx = idx[:np.min((len(edu), 1000))]
    edu = np.array(edu)
    plt.plot(edu[idx], log_edu[idx], '.c')
    px_edu = np.arange(-0.2, 5.3, .02)
    plt.plot(x_edu, np.max(px_edu) * (f_edu/np.max(f_edu)), 'k-')
    plt.plot(x_edu, np.max(px_edu) * (f2_edu/np.max(f2_edu)), 'k--')
    plt.plot([-200, 300, np.nan, -100, 500],
             [edu_s_min, edu_s_min, np.nan,
              edu_shift+np.log10(np.max((peak_loc-offset_edu, 1))),
              edu_shift+np.log10(np.max((peak_loc-offset_edu, 1)))], '-r')
    plt.plot([offset_edu, offset_edu], [0, 5], ':r')
    plt.plot([100, peak_trough, peak_trough],
             [edu_g1_max, edu_g1_max, 0], '--r')
    plt.ylim((px_edu[0], px_edu[-1]))
    plt.xlim([-200, np.max((peak_loc + 5 * peak_width, 500))])
    plt.xlabel('EdU intensity')
    plt.ylabel('log10 (EdU')
    return peak_loc, offset_edu, edu_g1_max, edu_s_min


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


def plot_edu_dna_scatter(dna, edu, offset_edu, x_dna=None, px_edu=None):
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)
    if px_edu is None:
        px_edu = np.arange(-0.2, 5.3, .02)
    log_dna = compute_log_dna(dna, x_dna)
    log_edu = compute_log_edu(edu, px_edu, offset_edu)
    xy = np.vstack([log_dna, log_edu])
    z = gaussian_kde(xy)(xy)
    plt.scatter(log_dna, log_edu, c=z, s=10)
    plt.xlabel('log10 (DNA)')
    plt.ylabel('log10 (EdU)')

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
    # Count the log intensity in each 2-D bin
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


def plot_2D_peaks(log_dna, x_dna, log_edu, px_edu, nsmooth=5):
    h = get_2d_histogram(log_dna, x_dna, log_edu, px_edu)
    g = smooth_1d(h, nsmooth)
    f = smooth_1d(g.T, nsmooth).T
    peak_candidates = iterate_2D_peak(h, x_dna, px_edu, nsmooth)
    plt.pcolor(x_dna, px_edu, f)
    plt.plot(peak_candidates[:, 0], peak_candidates[:, 1], 'ok',
             markersize=4, markerfacecolor='None')
    plt.xlabel('log (DNA)')
    plt.ylabel('log (EdU)')
    plt.xlim(quantile(x_dna, [0.25, 0.75]))
    plt.ylim(quantile(px_edu, [0.1, 0.75]))


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
            edu_peak_bool * np.sum(np.abs(repmat1 - repmat1.T < np.log10(2) * 0.75) * repmat2, axis=0).T +
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
def get_dna_peaks(log_dna, x_dna, log_edu, edu_shift,
                  min_edu, max_edu, phase_candidates):
    f_dna = get_kde(log_dna, x_dna)
    log_dna_low_edu = log_dna[(log_edu < min_edu+0.2*edu_shift) &
                              (log_edu < max_edu)]
    f_dna_low_edu = get_kde(log_dna_low_edu, x_dna)
    peak_amp, peak_loc, _ = findpeaks(f_dna_low_edu.tolist())
    peak_loc = peak_loc[peak_amp > np.max(peak_amp/10)]
    dna_peaks = x_dna[peak_loc[:3]]

    if len(dna_peaks) > 1:
        if phase_candidates[0, 0]:
            dna_peaks = dna_peaks[np.argmin(np.abs(
                dna_peaks - phase_candidates[0, 0]))]
        elif phase_candidates[1, 0]:
            dna_peaks = np.max(dna_peaks[
                dna_peaks > phase_candidates[1, 0]])
        else:
            dna_peaks = np.min(dna_peaks)

    if not np.any(dna_peaks):
        dna_peaks = np.nanmin(phase_candidates[:, 0] - np.log10(1.2))
    plt.plot(x_dna, f_dna)
    plt.plot(dna_peaks, .1, 'xk')
    plt.xlabel('log (DNA)')
    plt.ylabel('kernel density estimate')
    return dna_peaks
