import peakutils
import numpy as np
from scipy.stats import gaussian_kde
import scipy.io as sio
import matplotlib.pyplot as plt

def get_prominence_reference_level(signal, peak, peak_loc):
    
    left_range = signal[:peak_loc+1]
    if  peak != np.max(left_range):
        sorted_lr = sorted(left_range)
        sorted_peak_loc = sorted_lr.index(peak)
        left_crossing = sorted_lr[sorted_peak_loc + 1]
    else:
        left_crossing = left_range[0]
    # lc_index = signal.index(left_crossing)
    lc_indeces = np.where(np.array(signal) == left_crossing)[0]
    lc_index = [li for li in lc_indeces.tolist() if li < peak_loc][-1]
    # lc_index = min(lc_indeces, key=lambda x: abs(x-peak_loc))
    left_range_min = np.min(signal[lc_index:peak_loc+1])

    right_range = signal[peak_loc:]
    if peak != np.max(right_range):
        sorted_rr = sorted(right_range)
        sorted_peak_loc = sorted_rr.index(peak)
        right_crossing = sorted_rr[sorted_peak_loc + 1]
    else:
        right_crossing = right_range[-1]
    rc_indeces = np.where(np.array(signal) == right_crossing)[0]    
    rc_index = [ri for ri in rc_indeces if ri > peak_loc][0]
    # rc_index = signal.index(right_crossing)
    right_range_min = np.min(signal[peak_loc:rc_index+1])

    reference_level = np.max([left_range_min, right_range_min])
    reference_loc = signal.index(reference_level)

    return reference_loc, reference_level


def get_width_half_prominence(signal, peak, peak_loc):
    _, reference_level = get_prominence_reference_level(signal, peak, peak_loc)
    half_prominence = reference_level + 0.5 * (peak - reference_level)
    left_range = signal[:peak_loc+1]
    lhp = list(filter(lambda x: x<half_prominence, left_range))[-1]
    # lhp = min(left_range, key=lambda x:abs(x-half_prominence))
    left_width_indeces = np.where(np.array(signal) == (lhp))[0]
    left_width_ind = [li for li in left_width_indeces.tolist() if li < peak_loc][-1]
    right_range = signal[peak_loc:]
    rhp = list(filter(lambda x: x< half_prominence, right_range))[0]
    right_width_indeces = np.where(np.array(signal) == (rhp))[0]
    # rhp = min(right_range, key=lambda x:abs(x-half_prominence))
    right_width_ind = [ri for ri in right_width_indeces if ri > peak_loc][0]
    return left_width_ind, right_width_ind, half_prominence


def get_kde(x, x_grid, bandwidth=None):
    kde = gaussian_kde(x, bandwidth)
    return kde.evaluate(x_grid)


def findpeaks(signal):
    peak_loc = peakutils.peak.indexes(signal)
    peak_amp = [signal[loc] for loc in peak_loc]
    width = []
    for loc, value in zip(peak_loc, peak_amp):
        lw, lr, _ = get_width_half_prominence(signal, value, loc)
        width.append(lr-lw)
    return peak_amp, peak_loc, width    


# l = sio.loadmat('peaksig.mat')
# ps = l['PeakSig'][0].tolist()
# plt.plot(ps)
# peak_loc = peakutils.peak.indexes(ps)
# peak_values = [ps[ind] for ind in peak_loc]
# plt.scatter(peak_loc, peak_values)
# for loc, values in zip(peak_loc, peak_values):
#     ploc, pval = get_prominence_reference_level(ps, values, loc)
#     plt.plot([loc, loc], [pval, values])
#     lw, lr, hp = get_width_half_prominence(ps, values, loc)
#     plt.plot([lw, lr], [hp, hp])
    
