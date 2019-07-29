import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
import os
from matplotlib.collections import PolyCollection
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from cell_cycle_gating import dead_cell_filter as dcf



def plot_stacked_bar(data, ax=None,
                     data_cols=['G1', 'S', 'G2', 'S_dropout', 'other', 'M'],
                     condition_cols=['agent', 'cell_line', 'well'],
                     figname='Stack_bar_plot.png',
                     stacked=True,
                     title=None,
                     show_legend=True,
                     **kwargs):
    '''
    Make stacked barplot. Figure is saved to the current working folder.
    Parameters
    ========
    data: pd.DataFrame, report dataframe from the 'run_cell_cycle_gating' function
    ax: matplotlib.pyplot.Axes object, the axes to plot onto. 
    data_cols, list of str, columns to be plotted. By default it is ['G1', 'S', 'G2', 'S_dropout', 'other', 'M'].
    condition_cols: list or str, columns names to group the data, by default it is ['agent', 'cell_line', 'well'].
    figname: str, output figure name, if None, no figure will be saved. 
    stacked: bool, weather to plot stacked bars.
    '''
    ax = ax or plt.gca()
    if isinstance(condition_cols, str):
        data.set_index(condition_cols, inplace=True)
    else:
        data.set_index(data[condition_cols].apply(
            lambda x: '_'.join(x.values), axis=1).values, inplace=True)
    g = data[data_cols].plot(kind='bar', stacked=stacked, ax=ax, **kwargs)
    if show_legend:
        plt.legend(bbox_to_anchor=(1, 0.7))
    ax.tick_params(labelsize=16)
    ax.tick_params(axis='x', labelrotation=45)
    ax.axhline(0, color='black', linewidth=4)
    ax.set_ylabel('Fraction', fontsize=16)
    if title is not None:
        plt.title(title)
    if figname is None:
        return g
    else:
        plt.tight_layout()
        plt.savefig(figname)
        plt.close()


def batch_stacked_bar_plot(data,
                           data_cols=['G1', 'S', 'G2',
                                      'S_dropout', 'other', 'M'],
                           row_by='agent',
                           col_by='cell_line',
                           plot_x_col='concentration',
                           stacked=True,
                           savefig=False,
                           **kwargs):
    '''
    Make panels of stacked bars.
    Parameters
    ========
    data: pd.DataFrame, report dataframe from the 'run_cell_cycle_gating' function
    row_by, col_by: str, column names in the input data that separate subplots in rows and cols.
    plot_x_col: str, column name in the input data that set the x-axis.
    data_cols, list of str, columns to be plotted. By default it is ['G1', 'S', 'G2', 'S_dropout', 'other', 'M'].
    condition_cols: list or str, columns names to group the data, by default it is ['agent', 'cell_line', 'well'].
    '''

    data.sort_values([row_by, col_by, plot_x_col], inplace=True)
    nrows = data[row_by].unique().shape[0]
    ncols = data[col_by].unique().shape[0]
    fig, axes = plt.subplots(nrows, ncols, sharex=True, sharey=True, **kwargs)
    axes = axes.ravel()
    plot_idx = 0
    for i, row_value in enumerate(data[row_by].unique()):
        for j, col_value in enumerate(data[col_by].unique()):
            plot_data = data[(data[row_by] == row_value) &
                             (data[col_by] == col_value)]
            plot_stacked_bar(plot_data, ax=axes[plot_idx], data_cols=data_cols, condition_cols=plot_x_col,
                             figname=None, stacked=stacked, show_legend=False)
            axes[plot_idx].set_title('{}:{}'.format(
                str(row_value), str(col_value)))
            axes[plot_idx].set_xlabel(plot_x_col, fontsize=16)
            axes[plot_idx].legend([])
            plot_idx += 1

    if nrows >= ncols:
        legend_loc = (0.5, -0.25)
        plt.legend(loc='center', bbox_to_anchor=legend_loc,
                   fontsize='x-large', ncol=len(data_cols))
    else:
        legend_loc = (1.2, 0.5)
        plt.legend(loc='center', bbox_to_anchor=legend_loc,
                   fontsize='x-large')
    plt.ylabel('Fraction')
    if savefig:
        figname = '_'.join([condition, 'stacked_barplot.png'])
        plt.savefig(figname)
    else:
        plt.show()
    plt.close()
    return


def group_hist3d(dna_content, figsize=(8, 8), title=None, plot_num=0, figname=None):
    '''
    Make stacked histogram plots, each column is plotted as a series.
    '''
    #fig = plt.figure(figsize=figsize)
    #ax = plt.axes(projection='3d')
    #fig = plt.figure(figsize=(8, 4))
    grid_size=(2, 1)
    plot_num = plot_num % 2
    rel_pos = np.unravel_index(plot_num, grid_size)
    ax = plt.subplot2grid(grid_size, rel_pos, projection='3d')
    doses = list(dna_content.keys())
    verts = []
    for i, conc in enumerate(doses):
        hist_y, hist_x = np.histogram(dna_content[conc], 100, density=True)
        hist_x = hist_x[:-1]# + np.diff(hist_x[:2])
        verts.append(list(zip(hist_x, hist_y)))
        a = 0.1 * (i+1)
        #ax.bar3d(hist_x, i, 0, 0.1, 0.05, hist_y, color=(0.2, 0.7, 0.3, a), shade=False)#, alpha=0.2)
    poly = PolyCollection(verts, edgecolors='k')#, facecolors=[cc('r'), cc('g'), cc('b'), cc('y')])
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=np.arange(len(doses)), zdir='y')

    dna_content_values = []
    for k in dna_content.keys():
        dna_content_values += list(dna_content[k])
    
    ax.set_xlabel('log10(DNA content)')
    ax.set_xlim3d(4, 5)
    ax.set_ylabel('doses')
    ax.set_ylim3d(0, 9)
    #ax.set_zlabel('Z')
    ax.set_zlim3d(0, 5)
    ax.view_init(elev=67, azim=-45)
    ax.set_yticklabels(doses, fontsize=8)
    if title is not None:
        ax.set_title(title, fontsize=16)
    #if figname is not None:
    #    plt.savefig(figname)
    #else:
    #    return fig



 #def collate_wells(dfm, barcode):
 #    dfm = dfm[dfm.barcode == barcode].copy()
 #    conditions = list(zip(dfm.cell_line, dfm.agent))
 #    conditions = list(set(conditions))
 #    for condition in conditons:
def collate_wells(dfm, cell_line, agent, barcode=None):
     if barcode:
         dfm = dfm[dfm.barcode == barcode].copy()
     da = dfm[(dfm.cell_line == cell_line) & (dfm.agent == agent)][
             ['barcode', 'cell_line', 'agent', 'concentration', 'well']].copy()
     da = da.sort_values('concentration')
     dna_content = {}
     for conc in da.concentration.unique():
         wells = da[da.concentration == conc]['well'].tolist()
         barcodes = da[da.concentration == conc]['barcode'].tolist()
         wells = ["%s%s" % (w[0], int(w[1:])) for w in wells]
         dna_content[conc] = get_dna_content(barcodes, wells)
     return dna_content


def get_dna_content(barcodes, wells):
    dna_content = []
    for bc, well in zip(barcodes, wells):
        try:
            bc = [s for s in os.listdir('.') if bc in s][0]
        except IndexError:
            print("data folder with barcode %s does not exist" % bc)
        data_file = [s for s in os.listdir(bc) if 'result.%s[' % well in s][0]
        df = pd.read_table("%s/%s" % (bc, data_file))
        dna = np.array(df['Nuclei Selected - DNAcontent'].tolist())
        log_dna = dcf.compute_log_dna(dna)
        dna_content += list(log_dna)
    return np.array(dna_content)    



def plot(dfm, barcode=None, figname='test.pdf'):
    dfm = dfm[(dfm.cell_line.notnull()) & (dfm.agent.notnull())]
    if barcode:
         dfm = dfm[dfm.barcode == barcode].copy()     
    conditions = list(zip(dfm.cell_line, dfm.agent))
    conditions = list(set(conditions))
    nb_plots = len(conditions)
    pdf_pages = PdfPages(figname)
    nb_plots_per_page = 3

    for i, cond in enumerate(conditions):
        if i % nb_plots_per_page == 0:
            fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
        cell_line = cond[0]
        agent = cond[1]
        title =  "%s %s" % (cell_line, agent)
        dna_content = collate_wells(dfm, cell_line, agent)
        group_hist3d(dna_content, title=title, plot_num=i)
        
        if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
           plt.tight_layout()
           pdf_pages.savefig(fig)
           print("Completed analysis for %d out of %d conditions" %
                 (i+1, len(conditions)))
           plt.close('all')
    pdf_pages.close()       

             
         
                  
