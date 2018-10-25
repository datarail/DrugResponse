import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D


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


def group_hist3d(data, figsize=(8, 8), title=None, figname=None):
    '''
    Make stacked histogram plots, each column is plotted as a series.
    '''
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection='3d')
    for i in range(data.shape[1]):
        hist_y, hist_x = np.histogram(data.iloc[:, i], 100, density=True)
        hist_x = hist_x[:-1] + np.diff(hist_x[:2])
        ax.bar3d(hist_x, i, 0, 0.1, 0.05, hist_y, shade=False, alpha=0.5)
    if title is not None:
        ax.set_title(title, fontsize=16)
    ax.set_xlabel('DNA Content', fontsize=16)
    ax.set_zlabel('Density', fontsize=16)
    ax.set_yticks(np.arange(data.shape[1]))
    ax.set_yticklabels(data.columns.values, fontsize=16)
    ax.tick_params(axis='y', labelrotation=15)
    if figname is not None:
        plt.savefig(figname)
    else:
        return fig
