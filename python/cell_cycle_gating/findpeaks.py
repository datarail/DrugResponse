import peakutils
import numpy as np
from scipy.stats import gaussian_kde


def get_prominence_reference_level(signal, peak, peak_loc):
    """Returns the amplitude and location of the lower reference
    level of a peaks prominence.
    Note that prominence is the the length from the reference level
    upto the peak

    Paramters:
    ----------
    signal: 1D-array
    peak: float
       amplitude of peak whose prominece is to be computed
    peak_loc: float
       location on X-axis of peak whose prominence is to be computed

    Return:
    -------
    reference_loc : float
       location of X-axis of peak whose prominence is to be computed.
       Should equal peak_loc
    reference_level : float
       lower reference level of peak prominence.
    """
    left_range = signal[:peak_loc+1]
    if peak != np.max(left_range):
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
    """
    """
    _, reference_level = get_prominence_reference_level(signal, peak, peak_loc)
    half_prominence = reference_level + 0.5 * (peak - reference_level)
    left_range = signal[:peak_loc+1]
    lhp = list(filter(lambda x: x < half_prominence, left_range))[-1]
    # lhp = min(left_range, key=lambda x:abs(x-half_prominence))
    left_width_indeces = np.where(np.array(signal) == (lhp))[0]
    left_width_ind = [li for li in left_width_indeces.tolist()
                      if li < peak_loc][-1]
    right_range = signal[peak_loc:]
    rhp = list(filter(lambda x: x < half_prominence, right_range))[0]
    right_width_indeces = np.where(np.array(signal) == (rhp))[0]
    # rhp = min(right_range, key=lambda x:abs(x-half_prominence))
    right_width_ind = [ri for ri in right_width_indeces if ri > peak_loc][0]
    return left_width_ind, right_width_ind, half_prominence


def get_kde(x, x_grid, bandwidth=None):
    # https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
    if bandwidth:
        kde = gaussian_kde(x, bandwidth / x.std(ddof=1))
    else:
        kde = gaussian_kde(x)
    return kde.evaluate(x_grid)


def findpeaks(signal, npeaks=None):
    """Returns the amplitude , location and half-prominence width of peaks
    from the input signal

    Parameters:
    ----------
    signal : list
    npeaks : int
         number of peaks, locations and width returned sorted
         from highest to lowes amplitude peaks

    Returns:
    --------
    peak_amp : array of float
         amplitude of peaks
    peak_loc : array of float
         location of peaks
    width : array of float
        list of widths at half-prominence

    """
    peak_loc = peakutils.peak.indexes(signal, thres=0.25)
    peak_amp = [signal[loc] for loc in peak_loc]
    width = []
    for loc, value in zip(peak_loc, peak_amp):
        lw, lr, _ = get_width_half_prominence(signal, value, loc)
        width.append(lr-lw)
    if npeaks:
        sorted_peak_amps = sorted(peak_amp)[::-1]
        sorted_peak_indeces = [peak_amp.index(p) for p in sorted_peak_amps]
        sorted_loc = [peak_loc[i] for i in sorted_peak_indeces]
        sorted_width = [width[i] for i in sorted_peak_indeces]
        return (np.array(sorted_peak_amps[:npeaks]),
                np.array(sorted_loc[:npeaks]),
                np.array(sorted_width[:npeaks]))
    else:
        return np.array(peak_amp), np.array(peak_loc), np.array(width)
