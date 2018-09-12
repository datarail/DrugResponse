import pandas as pd
import os
import re
from cell_cycle_gating import ph3_filter as pf
from cell_cycle_gating import brdu_gating as bg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42



def get_gates_per_well(dfm, object_level_directory, object_level_files):
    """Runs BrdU and pH3 gating function scripts on a per well basis

    Parameters
    ----------
    dfm : pandas dataframe
        table containing well level metadata mapping wells to cell_line,
        agents, and concentration.
    object_level_directory : str
        path to folder containing object level .txt files
    object_level_files : list of str
        list of object level .txt files.

    Returns
    -------
    dfc : pandas dataframe
        metadata file appended with pH3 and BrdU gates.
    """
    # Merge agent and concentration columns in the case of
    # combination experiments.
    dfm = process_metadata_file(dfm)

    df_summary = pd.DataFrame()
    for file in object_level_files:
        gates = {}
        df = pd.read_table('%s/%s' % (object_level_directory, file))
        well = re.search('result.(.*?)\[', file).group(1)
        well = "%s%s" % (well[0], well[1:].zfill(2))
        gates['well'] = well

        brdu = np.array(df['Nuclei Selected - Nucleus A647 Mean'].tolist())
        ph3 = np.array(df['Nuclei Selected - Nucleus A488 Mean'].tolist())

        try:
            brdu_cutoff = bg.get_brdugates(brdu)
            _, ph3_cutoff, _ = pf.get_ph3_gates(ph3, cell_identity=None)
            ph3_cutoff = 10 ** ph3_cutoff
            gates['ph3_cutoff'] = ph3_cutoff
            gates['brdu_cutoff'] = brdu_cutoff
            df_summary = df_summary.append(gates, ignore_index=True)
        except IndexError:
            print(well)
            pass
    dfm.index = dfm['well']
    df_summary.index = df_summary['well']
    dfc = pd.concat([dfm, df_summary], axis=1)
    return dfc


def plot_scatter(dfc, dfm, object_level_directory,
                 object_level_files, filename='scatter.pdf'):
    """Plots BrdU and pH3 scatter and applies gating based on
        average gating of control wells

    Parameters
    ----------
    dfc : pandas dataframe
        table containing well level metadata and BrdU and pH3 gates per well.
    dfm : pandas dataframe
        table containing well level metadata mapping wells to cell_line,
        agents, and concentration.
    object_level_directory : str
        path to folder containing object level .txt files
    object_level_files : list of str
        list of object level .txt files.
    filename : Optional[str]
        Name of .pdf file in which the scatter plots for each well is saved.
        Default is 'scatter.pdf'.

    Returns
    -------
    dfc : pandas dataframe
        metadata file appended with pH3 and BrdU gates.
    """
    # Merge agent and concentration columns in the case of
    # combination experiments.
    dfm = process_metadata_file(dfm)

    # PDF setup
    pdf_pages = PdfPages(filename)
    nb_plots = len(object_level_files)
    nb_plots_per_page = 15
    # nb_pages = int(np.ceil(nb_plots / float(nb_plots_per_page)))
    plt.ioff()

    brdu_cutoff_mean = dfc.groupby(['agent', 'concentration'])[
        'brdu_cutoff'].mean().loc['DMSO'].values[0]
    ph3_cutoff_mean = dfc.groupby(['agent', 'concentration'])[
        'ph3_cutoff'].mean().loc['DMSO'].values[0]
    df_summary = pd.DataFrame()
    for i, file in enumerate(object_level_files):
        if i % nb_plots_per_page == 0:
            fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
        fractions = {}
        df = pd.read_table('%s/%s' % (object_level_directory, file))
        well = re.search('result.(.*?)\[', file).group(1)
        well = "%s%s" % (well[0], well[1:].zfill(2))
        agent = dfm[dfm.well == well].agent.values[0]
        concentration = dfm[dfm.well == well].concentration.values[0]
        if concentration == '0':
            conc = ''
        else:
            conc = "(%s um)" % concentration
            conc.replace('0.3333333333333', '0.33')
        title = "%s - %s %s" % (well, agent, conc)
        fractions['well'] = well

        brdu = np.array(df['Nuclei Selected - Nucleus A647 Mean'].tolist())
        # dna = np.array(df['Nuclei Selected - DNA Content'].tolist())
        ph3 = np.array(df['Nuclei Selected - Nucleus A488 Mean'].tolist())

        high_brdu = brdu[brdu > brdu_cutoff_mean]
        high_ph3 = ph3[ph3 > ph3_cutoff_mean]
        brduf = brdu[(brdu > brdu_cutoff_mean) & (ph3 > ph3_cutoff_mean)]
        ph3f = ph3[(brdu > brdu_cutoff_mean) & (ph3 > ph3_cutoff_mean)]
        fractions['total_cells'] = len(brdu)
        fractions['ph3_positive'] = len(high_ph3)
        fractions['brdu_positive'] = len(high_brdu)
        fractions['double_positive'] = len(brduf)
        fractions['brdu_cutoff_control_mean'] = brdu_cutoff_mean
        fractions['ph3_cutoff_control_mean'] = ph3_cutoff_mean
        fraction_s_phase = 100 * (len(brduf)/len(brdu))
        # fig = plt.figure()
        rel_pos = np.unravel_index(i % nb_plots_per_page, (5, 3))
        ax = plt.subplot2grid((5, 3), rel_pos)
        ax.scatter(brdu, ph3, color='black', s=10)
        ax.scatter(brduf, ph3f, color='red', s=10)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_title(title, fontsize=6)
        ax.set_xlabel('BrdU')
        ax.set_ylabel('pH3')
        # ylim = ax.get_ylim()
        # xlim = ax.get_xlim()
        ax.text(1.5e3, ph3_cutoff_mean+500, '%.2f%%' % fraction_s_phase,
                fontsize=10, color='black')
        ax.set_ylim(1e2, 2.5e4)
        ax.set_xlim(50, 5e3)
        ax.plot([brdu_cutoff_mean, brdu_cutoff_mean], [1e2, 2.5e4],
                'black', linewidth=1)
        ax.plot([50, 5e3], [ph3_cutoff_mean, ph3_cutoff_mean],
                'black', linewidth=1)
        df_summary = df_summary.append(fractions, ignore_index=True)

        if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
            plt.tight_layout()
            pdf_pages.savefig(fig)

    dfc.index = dfc['well']
    df_summary.index = df_summary['well']
    dfc = pd.concat([dfc, df_summary], axis=1)
    pdf_pages.close()
    return dfc


def plot_summary(dfm, object_level_folder, filename='scatter.pdf'):
    """Makes summary scatter plot of BrdU pH3 and returns summary
    data file of fraction of cells in each phase based on gating

    Parameters
    ----------
    dfm : pandas dataframe
        table containing well level metadata. Minimun columns required is
        well, cell_line, agent, and concentration.
    object_level_folder : str
        path to folder containing object level text files.
    filename : Optional[str]
         Name of .pdf file containing the summary plots per well.
         Default is scatter.pdf
    Returns
    -------
    dfs : pandas dataframe
        Summary of BrDu and pH3 based gating per well, and fraction of cells
        in each state.
    """
    dfm_ord = merge_metadata(dfm, object_level_folder)
    object_level_files = dfm_ord['object_level_file'].tolist()
    # object_level_files = [s for s in os.listdir(object_level_folder)
    #                       if 'Nuclei Selected[0].txt' in s]

    dfc = get_gates_per_well(dfm, object_level_folder, object_level_files)
    dfs = plot_scatter(dfc, dfm, object_level_folder,
                       object_level_files, filename)
    return dfs


def merge_metadata(dfm, obj):
    """Merges metadata with object level file and
    sorts files by cell line, drug and concetration.

    Parameters
    ----------
    dfm : pandas dataframe
        metadata file with wells mapped to treatment conditions
    obj : str
        path to folder containing object level data

    Returns
    -------
    dfmc : pandas dataframe
        mapping between wells, treatment, and filename
        containing object level data.
        Sorted by cell line, drug and concentration
    """
    # Merge agent and concentration columns in the case of
    # combination experiments.
    dfm = process_metadata_file(dfm)

    barcode = obj.split('[')[0]
    dfm = dfm[dfm.barcode == barcode].copy()
    dfm.index = dfm.well.tolist()
    wells = []
    object_level_files = [s for s in os.listdir(obj)
                          if 'Nuclei Selected[0].txt' in s]
    for file in object_level_files:
        well = re.search('result.(.*?)\[', file).group(1)
        well = "%s%s" % (well[0], well[1:].zfill(2))
        wells.append(well)
    dfmap = pd.DataFrame(list(zip(wells, object_level_files)),
                         columns=['well', 'object_level_file'])
    dfmap.index = dfmap.well.tolist()
    dfmc = pd.concat([dfm, dfmap], axis=1)
    dfmc = dfmc.dropna(subset=['object_level_file'])
    dfmc = dfmc.sort_values(['cell_line', 'agent', 'concentration'])
    return dfmc


def process_metadata_file(dfmeta):
    """Merge agent and concentration columns in the case of
       combination experiments.

    Parameters
    ----------
    dfmeta : pandas dataframe
       metadata file with wells mapped to treatment condition.
    Returns
    -------
    dfm : pandas dataframe
       metadata file with agent and concentration columns merged in the case of
    combination experiments.
    """
    dfm = dfmeta.copy()
    agent_cols = [a for a in dfm.columns.tolist() if 'agent' in a]
    if len(agent_cols) > 1:
        dfm[agent_cols] = dfm[agent_cols].replace(['', ' '], np.nan)
        dfm['agent'] = dfm[agent_cols].apply(
            lambda x: ','.join(x[x.notnull()]), axis=1)
        conc_cols = [a for a in dfm.columns.tolist() if 'concentration' in a]
        dfm[conc_cols] = dfm[conc_cols].astype(str)
        dfm[conc_cols] = dfm[conc_cols].replace(['0.0', '0'], np.nan)
        dfm['concentration'] = dfm[conc_cols].apply(
            lambda x: ','.join(x[x.notnull()]), axis=1)
    return dfm
