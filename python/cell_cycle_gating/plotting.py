import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from sklearn.preprocessing import robust_scale


def plot_stacked_bar(data, ax=None,
                     data_cols=['G1', 'S', 'G2', 'S_dropout', 'other', 'M'],
                     condition_cols=['agent', 'cell_line', 'well'],
                     figname='Stack_bar_plot.png',
                     stacked=True,
                     title=None,
                     show_legend=True,
                     **kwargs):
    """Make stacked barplot. Figure is saved to the current working folder.

    Parameters
    --------
    data: pd.DataFrame
        report dataframe from the 'run_cell_cycle_gating' function
    ax: matplotlib.pyplot.Axes
        the axes to plot onto.
    data_cols, list of str
        columns to be plotted. By default it is ['G1', 'S', 'G2', 'S_dropout', 'other', 'M'].
    condition_cols: list or str
        columns names to group the data, by default it is ['agent', 'cell_line', 'well'].
    figname: str
        output figure name, if None, no figure will be saved.
    stacked: bool
        weather to plot stacked bars.

    """
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
    """Make panels of stacked bars.

    Parameters
    --------
    data: pd.DataFrame
        report dataframe from the 'run_cell_cycle_gating' function.
    row_by, col_by: str
        column names in the input data that separate subplots in rows and cols.
    plot_x_col: str
        column name in the input data that set the x-axis.
    data_cols: list of str
        columns to be plotted. By default it is ['G1', 'S', 'G2', 'S_dropout', 'other', 'M'].
    condition_cols: list or str
        columns names to group the data, by default it is ['agent', 'cell_line', 'well'].
    """

    data.sort_values([row_by, col_by, plot_x_col], inplace=True)
    nrows = data[row_by].unique().shape[0]
    ncols = data[col_by].unique().shape[0]
    if ncols >= 2:
        total_plots = nrows * ncols
        ncols = 2
        nrows = int(np.ceil(total_plots / ncols))
    figsize = (ncols * 8, nrows * 6)
    fig, axes = plt.subplots(nrows, ncols, sharex=True,
                             sharey=True, figsize=figsize)
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
        legend_loc = (0.5, -0.5)
        plt.legend(loc='best', bbox_to_anchor=legend_loc,
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


def stacked_hist3d(data, title=None, figname=None, alpha=0.3, bins=50):
    """Make stacked histogram plots, each column is plotted as a series.

    Parameters
    --------
    title, figname: str
        name for title or figure, if None, no title/figure were included/saved.
    alpha: float
        transparency of each histogram.

    Returns
    --------
    fig:, matplotlib.pyplot.figure
        the figure object.

    """
    colors = sns.color_palette("husl", data.shape[1])
    if data.shape[1] * 2 < 8:
        figsize = (8, 8)
    else:
        figsize = (data.shape[1] * 2, data.shape[1] * 2)
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection='3d')
    # Add alpha values to the RGB tuples.
    colors = [x + (alpha,) for x in colors]
    for i in range(data.shape[1]):
        hist_data = data.iloc[:, i].dropna()
        if hist_data.max() > 10000:
            # hist_data = robust_scale(
            #     np.log2(hist_data + 1).values.reshape(-1, 1))
            hist_data = np.log10(hist_data + 1).values
            # print(hist_data[:5])
        hist_y, hist_x = np.histogram(hist_data, bins, density=True)
        hist_x = hist_x[:-1] + np.diff(hist_x[:2])
        # Use polygon3D to color the histogram, this is achieved by providing
        # all vertices of the polygon to be colored.
        verts = [(hist_x[j], i, hist_y[j]) for j in range(bins)] + \
            [(hist_x.max(), i, 0), (hist_x.min(), i, 0)
             ]  # Rememver to add the two vertices on the z=0 line.
        poly = Poly3DCollection([verts], color=colors[i])
        # poly.set_alpha(0.1) # This is not working currently.
        ax.add_collection3d(poly)
        ax.plot(hist_x, [i] * bins, hist_y, visible=False)
    if title is not None:
        ax.set_title(title, fontsize=16)
    ax.set_xlim(3, 7)
    ax.set_xlabel('DNA Content', fontsize=16)
    ax.set_zlabel('Density', fontsize=16)
    ax.set_yticks(np.arange(data.shape[1]))
    ax.set_yticklabels(data.columns.values, fontsize=16,
                       fontdict={'horizontalalignment': 'left'})
    ax.grid(alpha=0.3, visible=False)  # currently not working.
    # ax.tick_params(axis='y', labelrotation=-45)
    if figname is not None:
        plt.savefig(figname)
    else:
        return fig


def data_processing_stacked_hist3d(metadata, datapath, agents=[], doses=[], cells=[], pool_same_condition=True, plate_col='plate',
                                   well_col='well', drug_col='agent', dose_col='concentration', barcode_col='barcode', cell_col='cell_line'):
    """Data processing function for extracting data from the raw data folder. Where the root contains metadata of the experiment and each 
    sub folder has well level data. Drugs, doses and cells parameters are used to specify the conditions, or left blank so that the user
    can enter these later with a interactive prompt.

    Parameters
    --------
    metadata: pandas.DataFrame
        the metadata table in the root folder. Containing metadata specifying treament information.
    datapath: str or os.path object
        path to the root of raw data directory.
    agents, doses, cells: list of strings
        specify drug/dose/cell values to be extracted from raw data.
    pool_same_condition: bool
        determines if wells of the same condition are pooled or treated separately.
    (metadata)_col: str
        column names in the metadata table. If for any reason the metadata format are changed, these parameters can partially accomodate it.

    Returns
    --------
    processed_data: pandas.DataFrame
        a table that can be readily used in the stacked_hist3d function, each column is a series.
    """
    # Filtering based on metadata.
    processed_data = pd.DataFrame()
    for cond, cond_col in zip([agents, doses, cells], [drug_col, dose_col, cell_col]):
        if len(cond) == 0:
            all_conds = sorted(metadata[cond_col].astype(str).unique())
            # Get user inputs.
            conds = input(
                'Input {} value(s) or indices, separated by ",", from the list: '.format(cond_col) + ', '.join(all_conds))
            conds = conds.split(',')
            try:
                conds = [all_conds[int(x)] for x in conds]
            except:
                pass
        else:
            conds = cond
        metadata = metadata[metadata[cond_col].isin(conds)]
    # read DNA content data iteratively based on the filtered metadata.
    target_dirs = metadata[barcode_col].unique()
    # iterate through plates
    for _dir in target_dirs:
        dir_metadata = metadata[metadata[barcode_col] == _dir]
        target_wells = dir_metadata[well_col].unique()
        dir_name = [x for x in os.listdir(datapath) if (
            (_dir in x) & ('.' not in x))][0]
        data_fn_list = os.listdir(os.path.join(datapath, dir_name))
        # iterate through wells.
        for _well in target_wells:
            agent_name = dir_metadata.set_index(well_col)[drug_col][_well]
            dose_name = dir_metadata.set_index(well_col)[dose_col][_well]
            cell_name = dir_metadata.set_index(well_col)[cell_col][_well]
            # File name and well name is not the same, temporary fix.
            if _well[1] == '0':
                _well = _well.replace('0', '')
            anchor_string = '.' + _well

            fn = [x for x in data_fn_list if anchor_string in x][0]
            target_data = pd.read_table(os.path.join(datapath, dir_name, fn))
            dna_content = target_data[
                [x for x in target_data.columns if 'DNAcontent' in x]]
            # Pool wells of same conditions.
            if pool_same_condition:
                dna_content.columns = [
                    '_'.join([agent_name, dose_name, cell_name])]
            else:
                dna_content.columns = [
                    '_'.join([agent_name, dose_name, cell_name, _well])]
            # Concatenate data column-wise, if column already exists, append
            # data instead..
            if dna_content.columns[0] in processed_data.columns:
                processed_data = processed_data.append(
                    dna_content, ignore_index=True)
            else:
                processed_data = pd.concat(
                    [processed_data, dna_content], axis=1)
    return processed_data
