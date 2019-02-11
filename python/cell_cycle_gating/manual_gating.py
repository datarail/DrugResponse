import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.widgets import Slider, Button
import numpy as np
import matplotlib.pyplot as plt
from cell_cycle_gating import cellcycle_phases as cc
from cell_cycle_gating import dead_cell_filter as dcf
from cell_cycle_gating import ph3_filter as pf
import pandas as pd
import re


def update_gating(dna, edu, px_edu=None, x_dna=None):
    if px_edu is None:
        px_edu = np.arange(-0.2, 5.3, .02)
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)

    log_dna = cc.compute_log_dna(dna, x_dna)

    edu_shift, offset_edu, edu_g1_max, edu_s_min = cc.get_edu_gates(edu, px_edu)
    
    log_edu = cc.compute_log_edu(edu, px_edu, offset_edu)

    dna_gates = [4.5, 5,  5.6, 6]
    edu_gates = [2, 3.5]

    xy = np.vstack([log_dna, log_edu])
    z = cc.gaussian_kde(xy)(xy)

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.5)
    plt.scatter(log_dna, log_edu, c=z, s=2, rasterized=True)
    l, = plt.plot([dna_gates[i] for i in [0, 0, 3, 3, 0, 0, 1, 1, 0, 2, 2, 3]],
                  [-1, edu_gates[1], edu_gates[1], -1,
                   np.nan, edu_gates[0], edu_gates[0],
                   -1, np.nan, -1, edu_gates[0], edu_gates[0]],
                  '--',  color='red')
    plt.xlabel('log10 (DNA content)')
    plt.ylabel('log10 (EdU)')
    plt.xlim((4, 6))
    plt.ylim((0, 4))
    #plt.xlim((log_dna.min(), log_dna.max()))
    #plt.ylim((log_edu.min(), log_edu.max()))
    axcolor = 'lightgoldenrodyellow'
    edu_lb = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    edu_ub = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    sedu_lb = Slider(edu_lb, 'EdU lower bound', 0.1, 10.0,
                     valinit=edu_gates[0])
    sedu_ub = Slider(edu_ub, 'Edu upper bound', 0.1, 10.0,
                     valinit=edu_gates[1])

    g1_left = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
    g1_right = plt.axes([0.25, 0.25, 0.65, 0.03], facecolor=axcolor)
    g2_left = plt.axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)
    g2_right = plt.axes([0.25, 0.35, 0.65, 0.03], facecolor=axcolor)

    sg1_left = Slider(g1_left, 'G1 left bound', 0.1, 10.0,
                     valinit=dna_gates[0])
    sg1_right = Slider(g1_right, 'G1 right bound', 0.1, 10.0,
                     valinit=dna_gates[1])
    sg2_left = Slider(g2_left, 'G2 left bound', 0.1, 10.0,
                     valinit=dna_gates[2])
    sg2_right = Slider(g2_right, 'G2 right bound', 0.1, 10.0,
                     valinit=dna_gates[3])


    def update(val):
        nedu_lb = sedu_lb.val
        nedu_ub = sedu_ub.val
        l.set_ydata([-1, nedu_ub, nedu_ub, -1,
                     np.nan, nedu_lb, nedu_lb,
                     -1, np.nan, -1, nedu_lb, nedu_lb])

        ndna_gates = [sg1_left.val, sg1_right.val, sg2_left.val, sg2_right.val]
        l.set_xdata([ndna_gates[i] for i in [0, 0, 3, 3, 0, 0, 1, 1, 0, 2, 2, 3]])
        fig.canvas.draw_idle()
    sedu_lb.on_changed(update)
    sedu_ub.on_changed(update)
    sg1_left.on_changed(update)
    sg1_right.on_changed(update)
    sg2_left.on_changed(update)
    sg2_right.on_changed(update)

    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


    def reset(event):
        sedu_lb.reset()
        sedu_ub.reset()
        sg1_left.reset()
        sg1_right.reset()
        sg2_left.reset()
        sg2_right.reset()
    button.on_clicked(reset)
    #global close_event
    #close_event = 0
    # def reset2(event):
    #     global close_event
    #     #fig.canvas.mpl_disconnect(cid)
    #     close_event += 1
    #     print(close_event)
    #     plt.close()

    #print('after', close_event)    
    #button.on_clicked(reset2)
    #if close_event:
    #    plt.pause(1)
    #else:
    plt.pause(30)

    #def onclick(event):
    #    global offset
    #    offset = event.ydata
    #    fig.canvas.mpl_disconnect(cid)
    #    plt.close()
    #    return     
        
    

    #cid = fig.canvas.mpl_connect('button_press_event', onclick)
    #plt.pause(30)
     
    #plt.close()

    
    ndna_gates = [sg1_left.val, sg1_right.val, sg2_left.val, sg2_right.val]
    nedu_gates = [sedu_lb.val, sedu_ub.val]

    #fractions = reevaluate_phases(log_dna, ndna_gates, log_edu, nedu_gates)

    return ndna_gates, nedu_gates


def reevaluate_phases(log_dna, dna_gates, log_edu, edu_gates):
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
    return fractions, cell_id


def manual_gating(dfs, obj_file, ndict, ldr_channel=True, ph3_channel=True,
                  x_dna=None, px_edu=None):
    df = pd.read_table(obj_file)
    df = map_channel_names(df, ndict)
    well = re.search('result.(.*?)\[', obj_file).group(1)
    well = "%s%s" % (well[0], well[1:].zfill(2))

    #edu = np.array(df['Nuclei Selected - EdUINT'].tolist())
    #dna = np.array(df['Nuclei Selected - DNAcontent'].tolist())
    #ldr = np.array(df['Nuclei Selected - LDRTXT SER Spot 8 px'].tolist())
    edu = np.array(df['edu'].tolist())
    dna = np.array(df['dna'].tolist())

    edu_notnan = ~np.isnan(edu)
    edu = edu[edu_notnan]
    dna = dna[edu_notnan]

    if ldr_channel:
        ldr = np.array(df['ldr'].tolist())
        ldr = ldr[edu_notnan]     

    if ph3_channel:
        #ph3 = np.array(df['Nuclei Selected - pH3INT'].tolist())
        ph3 = np.array(df['ph3'].tolist())
        ph3 = ph3[edu_notnan]

    if px_edu is None:
        px_edu = np.arange(-0.2, 5.3, .02)
    if x_dna is None:
        x_dna = np.arange(2.5, 8, 0.02)  
    log_dna = cc.compute_log_dna(dna, x_dna)
    edu_shift, offset_edu, edu_g1_max, edu_s_min = cc.get_edu_gates(edu, px_edu) 
    log_edu = cc.compute_log_edu(edu, px_edu, offset_edu)    

    ndna_gates, nedu_gates = update_gating(dna, edu)    
    fractions, cell_id = reevaluate_phases(log_dna, ndna_gates, log_edu, nedu_gates)    

    if ldr_channel:
        ldr_gates = dcf.get_ldrgates(ldr)
        odna_gates = dcf.get_dna_gating(dna, ldr, ldr_gates)
        a, d, _ = dcf.live_dead(ldr, ldr_gates, dna, odna_gates)    

    if ph3_channel:
        f_ph3, ph3_cutoff, ph3_lims = pf.get_ph3_gates(ph3, cell_id)
        log_ph3 = pf.compute_log_ph3(ph3)
        fractions = pf.evaluate_Mphase(log_ph3, ph3_cutoff, cell_id)

    if ldr_channel:
        fractions['cell_count'] = a
        fractions['cell_count__dead'] = d
    fractions['well'] = well
    fractions['cell_count__total'] = len(dna)

    if 'corpse_count' in dfs.columns.tolist():
        fractions['corpse_count'] = dfs[dfs.well == well]['corpse_count'].values[0]

    #dfs2 = dfs[dfs.well != well].copy()
    dnew = pd.DataFrame(fractions, index=[well])
    dfs2 = dfs.copy()
    dfs2.index = dfs2['well']
    #dfs2 = dfs2.append(fractions, ignore_index=True)
    dfs2.update(dnew)
    return dfs2


def map_channel_names(df, ndict):
    df = df.rename(columns=ndict)
    return df
