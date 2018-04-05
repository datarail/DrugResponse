import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re
from matplotlib.backends.backend_pdf import PdfPages
from cell_cycle_gating import dead_cell_filter as dcf
from cell_cycle_gating import cellcycle_phases as cc
from cell_cycle_gating import ph3_filter as pf

# object_level_directory = '180202_11II_062[20634]'


def run(object_level_directory):
    plt.ioff()
    object_level_files = [s for s in os.listdir(object_level_directory)
                          if 'Nuclei Selected[0].txt' in s]

    df_summary = pd.DataFrame()
    nb_plots = len(object_level_files)
    pdf_pages = PdfPages('summary_%s.pdf' % object_level_directory)
    nb_plots_per_page = 10
    # nb_pages = int(np.ceil(nb_plots / float(nb_plots_per_page)))
    for i, file in enumerate(object_level_files[:nb_plots]):
        if i % nb_plots_per_page == 0:
            fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
        df = pd.read_table('%s/%s' % (object_level_directory, file))
        well = re.search('result.(.*?)\[', file).group(1)

        edu = np.array(df['Nuclei Selected - EdUINT'].tolist())
        dna = np.array(df['Nuclei Selected - DNAcontent'].tolist())
        ldr = np.array(df['Nuclei Selected - LDRTXT SER Spot 8 px'].tolist())
        ph3 = np.array(df['Nuclei Selected - pH3INT'].tolist())

        edu_notnan = ~np.isnan(edu)
        edu = edu[edu_notnan]
        dna = dna[edu_notnan]
        ldr = ldr[edu_notnan]
        ph3 = ph3[edu_notnan]

        try:
            # Get live dead
            ldr_gates = dcf.get_ldrgates(ldr)
            dna_gates = dcf.get_dna_gating(dna, ldr, ldr_gates)
            a, d, _ = dcf.live_dead(ldr, ldr_gates, dna, dna_gates)
            # Get phases based on DNA and EdU
            fractions, cell_identity = cc.plot_summary(dna, edu, fig,
                                                       well=well,
                                                       plot='scatter',
                                                       plot_num=i)

            # Revise phases based on pH3
            f_ph3, ph3_cutoff, ph3_lims = pf.get_ph3_gates(ph3, cell_identity)
            log_ph3 = pf.compute_log_ph3(ph3)
            fractions = pf.evaluate_Mphase(log_ph3, ph3_cutoff, cell_identity)

            fractions['well'] = well
            fractions['alive'] = a
            fractions['dead'] = d
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
        if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
            plt.tight_layout()
            pdf_pages.savefig(fig)
    pdf_pages.close()
    df_summary = df_summary[['well', 'alive', 'dead',
                             'G1', 'S', 'G2', 'M',
                             'S_dropout', 'other']]
    df_summary.to_csv('summary_%s.csv' % object_level_directory, index=False)
    return df_summary
