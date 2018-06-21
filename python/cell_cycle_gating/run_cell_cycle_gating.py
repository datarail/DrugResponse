import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re
from matplotlib.backends.backend_pdf import PdfPages
from cell_cycle_gating import dead_cell_filter as dcf
from cell_cycle_gating import cellcycle_phases as cc
from cell_cycle_gating import ph3_filter as pf


def run(object_level_directory, dfm=None, ph3_channel=True):
    """ Executes cell cycle gating on all wells for which object level
    data is available in object_level_directory. Plots and saves summary pdf
    of DNA v EdU distribution with automated gatings. A dataframe summarizing
    results for each well is also saved.

    Parameters
    ----------
    object_level_directory : str
        name of folder containing object level data for a single plate
    dfm : Optional[pandas dataframe]
        metadata table of experimental design. Default is None.
    ph3_channel : Optional[bool]
        True if data on pH3 intensity is in object level data.
    Returns
    -------
    df : pandas dataframe
       well level summary of number of live/dead cells and fraction of
       cells in each phase of the cell cycle
    """
    plt.ioff()
    if dfm is not None:
        dfm_ord = merge_metadata(dfm, object_level_directory)
        object_level_files = dfm_ord['object_level_file'].tolist()
    else:
        object_level_files = [s for s in os.listdir(object_level_directory)
                              if 'Nuclei Selected[0].txt' in s]

    df_summary = pd.DataFrame()
    nb_plots = len(object_level_files)
    pdf_pages = PdfPages('summary_%s.pdf' % object_level_directory)
    nb_plots_per_page = 10
    # nb_pages = int(np.ceil(nb_plots / float(nb_plots_per_page)))
    for i, file in enumerate(object_level_files):
        if i % nb_plots_per_page == 0:
            fig = plt.figure(figsize=(8.27, 11.69), dpi=10)
        df = pd.read_table('%s/%s' % (object_level_directory, file))
        well = re.search('result.(.*?)\[', file).group(1)
        well = "%s%s" % (well[0], well[1:].zfill(2))

        edu = np.array(df['Nuclei Selected - EdUINT'].tolist())
        dna = np.array(df['Nuclei Selected - DNAcontent'].tolist())
        ldr = np.array(df['Nuclei Selected - LDRTXT SER Spot 8 px'].tolist())
        
        edu_notnan = ~np.isnan(edu)
        edu = edu[edu_notnan]
        dna = dna[edu_notnan]
        ldr = ldr[edu_notnan]

        if ph3_channel:
            ph3 = np.array(df['Nuclei Selected - pH3INT'].tolist())
            ph3 = ph3[edu_notnan]

        try:
            # Get live dead
            ldr_gates = dcf.get_ldrgates(ldr)
            dna_gates = dcf.get_dna_gating(dna, ldr, ldr_gates)
            a, d, _ = dcf.live_dead(ldr, ldr_gates, dna, dna_gates)
            # Get phases based on DNA and EdU
            if dfm is not None:
                # dfm_ord.index = dfm_ord.well
                cell_line = dfm_ord.loc[well, 'cell_line']
                agent = dfm_ord.loc[well, 'agent']
                conc = dfm_ord.loc[well, 'concentration']
                title = "%s %s %s (%s um)" % (well, cell_line, agent, conc)
            else:
                title = well
            fractions, cell_identity = cc.plot_summary(dna, edu, fig,
                                                       title=title,
                                                       plot='scatter',
                                                       plot_num=i)

            # Revise phases based on pH3
            if ph3_channel:
                f_ph3, ph3_cutoff, ph3_lims = pf.get_ph3_gates(ph3, cell_identity)
                log_ph3 = pf.compute_log_ph3(ph3)
                fractions = pf.evaluate_Mphase(log_ph3, ph3_cutoff, cell_identity)

            fractions['well'] = well
            fractions['cell_count__total'] = len(dna)
            fractions['cell_count'] = a
            fractions['cell_count__dead'] = d
            df_summary = df_summary.append(fractions, ignore_index=True)
        except ValueError:
            print(well, ' ValueError')
            pass
        except TypeError:
            print(well, 'TypeError')
            pass
        except IndexError:
            print(well, 'IndexError')
            pass
    # your code that will (maybe) throw
        except np.linalg.LinAlgError as e:
            if 'Singular matrix' in str(e):
                print(well, 'Singular matrix error')
                pass
        if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
            plt.tight_layout()
            pdf_pages.savefig(fig)
            print("Completed analysis for %d out of %d wells" %
                  (i+1, len(object_level_files)))
            plt.close('all')
    pdf_pages.close()
    if ph3_channel:
        summary_cols = ['well', 'cell_count__total',
                        'cell_count', 'cell_count__dead',
                        'G1', 'S', 'G2', 'M',
                        'S_dropout', 'other']
    else:
        summary_cols = ['well', 'cell_count__total',
                        'cell_count', 'cell_count__dead',
                        'G1', 'S', 'G2',
                        'S_dropout', 'other']
    df_summary = df_summary[summary_cols]
    # Merge summary table with metadata if provided
    if dfm is not None:
        df_summary.index = df_summary['well'].tolist()
        df_summary = pd.concat([dfm_ord, df_summary], axis=1)
        df_summary = df_summary.loc[:, ~df_summary.columns.duplicated()]
    # Merge summary table with corpse count if available
    dfc = get_corpse_count(object_level_directory)
    if dfc is not None:
        df_summary.index = df_summary['well'].tolist()
        df_summary = pd.concat([df_summary, dfc], axis=1)
    df_summary.to_csv('summary_%s.csv' % object_level_directory, index=False)
    return df_summary


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


def get_corpse_count(obj):
    """Return corpse count dataframe

    Parameters
    ----------
    obj : str
        path to folder containing object level data

    Returns
    -------
    dfc : pandas dataframe
        datatable mapping well to corpse counts
    """
    corpse_file = [f for f in os.listdir(obj)
                   if 'CorpseCountOnly' in f and '.txt' in f]
    if corpse_file:
        df_corpse = pd.read_table("%s/%s" %
                                  (obj, corpse_file[0]))
        wells = ["%s%s" % (w[0], w[1:].zfill(2))
                 for w in df_corpse.WellName.tolist()]
        df_corpse.index = wells
        nd = {'Corpses - Number of Objects': 'corpse_count'}
        df_corpse = df_corpse.rename(columns=nd)
        df_corpse = df_corpse['corpse_count'].copy()
        dfc = pd.DataFrame(df_corpse)
        return dfc
    else:
        return None
