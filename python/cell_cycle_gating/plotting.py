import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_stacked_bar(data, ax=None,
    data_cols = ['G1', 'S', 'G2', 'S_dropout','other', 'M'],
    condition_cols = ['agent','cell_line','well'],
    figname = 'Stack_bar_plot.png',
    stacked = True,
    title = None,
    show_legend = True,
    **kwargs):
    '''
    Make stacked barplot. Figure is saved to the current working folder.
    Parameters
    ========
    Todo
    '''
    ax = ax or plt.gca()
    if isinstance(condition_cols,str):
        data.set_index(condition_cols, inplace=True)
    else:
        data.set_index(data[condition_cols].apply(lambda x: '_'.join(x.values),axis = 1).values, inplace=True)
    g = data[data_cols].plot(kind='bar', stacked=stacked, ax=ax, **kwargs)
    if show_legend:
        plt.legend(bbox_to_anchor = (1,0.7))
    ax.tick_params(labelsize=16)
    ax.tick_params(axis = 'x',labelrotation=45)
    plt.ylabel('Fraction')
    if title is not None:
        plt.title(title)
    if figname is None:
        return g
    else:
        plt.tight_layout()
        plt.savefig(figname)
        plt.close()

def batch_stacked_bar_plot(data,
    data_cols = ['G1', 'S', 'G2', 'S_dropout','other', 'M'],
    row_by = 'agent',
    col_by = 'cell_line',
    plot_x_col = 'concentration',
    stacked = True,
    savefig = False,
    **kwargs):
    data.sort_values([row_by,col_by,plot_x_col],inplace=True)
    nrows = data[row_by].unique().shape[0]
    ncols = data[col_by].unique().shape[0]
    fig, axes = plt.subplots(nrows,ncols,sharex=True, sharey=True,**kwargs)
    axes = axes.ravel()
    plot_idx = 0
    for i,row_value in enumerate(data[row_by].unique()):
        for j,col_value in enumerate(data[col_by].unique()):
            plot_data = data[(data[row_by]==row_value) & (data[col_by]==col_value)]
            plot_stacked_bar(plot_data, ax=axes[plot_idx], data_cols=data_cols, condition_cols=plot_x_col, figname=None,stacked=stacked, show_legend = False)
            axes[plot_idx].set_title('{}:{}'.format(str(row_value),str(col_value)))
            axes[plot_idx].legend([])
            plot_idx+=1
    if nrows>=ncols:
        legend_loc = (1.1,-0.5)
        plt.legend(bbox_to_anchor = legend_loc,ncol=len(data_cols))
    else:
        legend_loc = (1,0.7)
        plt.legend(bbox_to_anchor = legend_loc)
    plt.ylabel('Fraction')
    if savefig:
        figname = '_'.join([condition,'stacked_barplot.png'])
        plt.savefig(figname)
    else:
        plt.show()
    plt.close()
    return