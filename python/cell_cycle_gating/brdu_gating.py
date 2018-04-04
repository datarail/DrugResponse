from cell_cycle_gating import findpeaks
import numpy as np
import matplotlib.pyplot as plt


def get_brdugates(brdu, x_brdu=None, plotting=False):
    if x_brdu is None:
        mx = np.max(brdu.tolist())+0.01
        x_brdu = np.arange(-0.01, mx, 1)
    f_brdu = findpeaks.get_kde(brdu, x_brdu)  # brdu should be an array
    peak_amp, peak_loc, peak_width = findpeaks.findpeaks(
        f_brdu.tolist(), npeaks=1)

    # choose BRDU cutoff based on half-proximal width and
    # right trough of peak
    width_2p5 = int((peak_loc + 2.5 * peak_width[0])[0])
    width_5 = int((peak_loc + 5 * peak_width[0])[0])

    # Find location of minimun on the right
    f_neg = [-x for x in f_brdu[width_2p5:width_5]]
    _, trough_loc, _ = findpeaks.findpeaks(f_neg, npeaks=1)
    if np.any(trough_loc):
        trough_loc = trough_loc[0] + peak_loc[0] - 1
    else:
        trough_loc = width_2p5
    brdu_cutoff = x_brdu[trough_loc]
    if plotting:
        plt.plot(x_brdu, f_brdu)
        plt.plot([brdu_cutoff, brdu_cutoff],
                 [0, 0.5 * peak_amp])
    return brdu_cutoff
