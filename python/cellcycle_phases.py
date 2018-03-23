from findpeaks import findpeaks, get_kde
import numpy as np
from scipy.stats.mstats import mquantiles as quantile
import matplotlib.pyplot as plt
import smooth
from scipy.stats import gaussian_kde
import math
import pandas as pd

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
    # Note: Bandwidth = 100 reproduced MATLAB outpu
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
        peak_trough = peak_loc + 3*peak_width

    # Edu offset
    # ** Not entirely clear to me yet
    offset_edu = np.max((peak_loc-1.5 * peak_width, 1))

    log_edu = compute_log_edu(edu, px_edu, offset_edu)

    # EdU max for G1 gating 
    edu_g1_max = np.max((
        # Expected EdU max for G1 (option 1)
        np.log10(peak_trough - offset_edu),
        quantile(log_edu, 0.2) + 0.1  # Expected EdU max  for G1 (option 2)
    ))

    # Edu  min for S phase gating
    edu_s_min = np.max((
        # Expected EdU min for S phase (option 1)
        np.log10(peak_loc + 2 * peak_width - offset_edu),
        edu_g1_max - 0.1  # Expected Edu min for S phase (option 2)
        ))

    # Expected differene between G1/G2 and S
    edu_shift = np.max((
        # Expected EdU min for S phase (option 1) - G1 peak location
        (np.log10(peak_loc + 2 * peak_width - offset_edu) -
         np.log10(np.max((peak_loc - offset_edu, 1)))),
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


def plot_edu_gating_estimate(edu, px_edu=None):
    if not px_edu:
        px_edu = np.arange(-0.2, 5.3, .02)
    edu_g1_max, edu_s_min, edu_shift, offset_edu = get_edu_gates(edu)


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
