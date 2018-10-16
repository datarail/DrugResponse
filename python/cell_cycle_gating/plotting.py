import pandas as pd
import matplotlib.pyplot as plt

def plot_stacked_bar(data,
    data_cols = ['G1', 'S', 'G2', 'S_dropout','other', 'M'],
    metadata_cols = ['agent','cell_line','well'],
    figname = 'Stack_bar_plot.png',
    stacked = True,
    title = None):
    '''
    Make stacked barplot. Figure is saved to the current working folder.
    Parameters
    ========
    Todo
    '''
    if isinstance(metadata_cols,str):
        data.set_index(metadata_cols, inplace=True)
    else:
        data.set_index(data[metadata_cols].apply(lambda x: '_'.join(x.values),axis = 1).values, inplace=True)
    g = data[data_cols].plot(kind='bar', stacked=stacked)
    plt.legend(bbox_to_anchor = (1,0.7))
    plt.xticks(rotation=45)
    plt.ylabel('Fraction')
    if title is not None:
        plt.title(title)
    if figname is None:
        plt.show()
    else:
        plt.tight_layout()
        plt.savefig(figname)
        plt.close()

def batch_stacked_bar_plot(data,
    data_cols = ['G1', 'S', 'G2', 'S_dropout','other', 'M'],
    condition_col = ['agent','cell_line','well'],
    plot_x_col = 'concentration',
    stacked = True):
    data['condition'] = data[condition_col].apply(lambda x: '_'.join(x.values), axis = 1)
    all_conditions = data.condition.unique()
    for condition in all_conditions:
        condition_data = data[data.condition==condition].iloc[:,:-1]
        figname = '_'.join([condition,'stacked_barplot.png'])
        condition_data.sort_values(plot_x_col,inplace=True)
        plot_stacked_bar(condition_data,data_cols=data_cols, metadata_cols=plot_x_col, figname=figname,stacked=stacked, title=condition)
