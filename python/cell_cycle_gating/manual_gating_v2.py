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
from ipywidgets import interactive, interact, fixed, interact_manual
import ipywidgets as widgets
from IPython.display import display
import os


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


def update_gating(dfs, obj, well, ndict,
                  ldr_channel=True, ph3_channel=True,
                  x_dna=None, px_edu=None):
    obj_file = get_obj_file(obj, well)
    path2file = "%s/%s" % (obj, obj_file)
    df = pd.read_table(path2file)
    df = map_channel_names(df, ndict)
    well = re.search('result.(.*?)\[', obj_file).group(1)
    well = "%s%s" % (well[0], well[1:].zfill(2))

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

    y=interactive(gating,
                  log_dna=fixed(log_dna), 
                  g1_left= np.median(log_dna)-.25,
                  g1_right = np.median(log_dna) - 0.15,
                  g2_left= np.median(log_dna) + 0.15,
                  g2_right = np.median(log_dna) +.25, 
                  log_edu = fixed(log_edu),
                  edu_lower = np.median(log_edu),
                  edu_upper = np.median(log_edu) + 0.5)
    for i, child in enumerate(y.children):
        child.step = 0.05
        if i <=3:
            child.min = log_dna.min()
            child.max = log_dna.max()
        else:
            child.min = log_edu.min()
            child.max = log_edu.max()
    return y


def gating(log_dna, log_edu,
            g1_left, g1_right, g2_left, g2_right, edu_lower, edu_upper):
    dna_gates = [g1_left, g1_right, g2_left, g2_right]
    edu_gates = [edu_lower, edu_upper]    
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
    plt.xlim((log_dna.min(), log_dna.max()))
    plt.ylim((log_edu.min(), log_edu.max()))
    
    #plt.pie(fractions.values(), labels=fractions.keys())
    plt.show()
    

def apply_gating(y, dfs, obj, well, ndict,
                  ldr_channel=True, ph3_channel=True,
                  x_dna=None, px_edu=None):
    obj_file = get_obj_file(obj, well)
    path2file = "%s/%s" % (obj, obj_file)
    df = pd.read_table(path2file)
    df = map_channel_names(df, ndict)
    well = re.search('result.(.*?)\[', obj_file).group(1)
    well = "%s%s" % (well[0], well[1:].zfill(2))

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

    ndna_gates = [y.kwargs['g1_left'], y.kwargs['g1_right'], y.kwargs['g2_left'], y.kwargs['g2_right']]
    nedu_gates = [y.kwargs['edu_lower'], y.kwargs['edu_upper']]
     
    fractions, cell_id = reevaluate_phases(log_dna, ndna_gates, log_edu, nedu_gates)    

    if ph3_channel:
        f_ph3, ph3_cutoff, ph3_lims = pf.get_ph3_gates(ph3, cell_id)
        log_ph3 = pf.compute_log_ph3(ph3)
        fractions = pf.evaluate_Mphase(log_ph3, ph3_cutoff, cell_id)

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


def get_obj_file(obj, well):
    if well[1] == '0':
        well_name = "%s%s" % (well[0], well[2:])
    else:
        well_name = well
    files = [f for f in os.listdir(obj)
             if f.endswith('Nuclei Selected[0].txt')]
    obj_file = [file for file in files if well_name ==
                re.search('result.(.*?)\[', file).group(1)]
    return obj_file[0]
