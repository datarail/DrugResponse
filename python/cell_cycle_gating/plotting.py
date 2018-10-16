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

def median_std_data(data, 
    drug_col='agent', 
    time_col='time', 
    cell_col='cell_line', 
    dose_col='concentration', 
    gr_col='GRvalue'):
    median_data = data.groupby([drug_col,time_col,cell_col,dose_col])[gr_col].median()
    std_data = data.groupby([drug_col,time_col,cell_col,dose_col])[gr_col].std()
    processed_data = pd.concat([median_data,std_data],axis = 1)
    processed_data.columns = ['median_GR','std_GR']
    processed_data.sort_values('median_GR',ascending=False,inplace=True)
    return processed_data

def waterfall_plot_process_data(separated_data, 
    condense_by = [''], 
    drug_col='agent', 
    time_col='time', 
    cell_col='cell_line', 
    dose_col='concentration', 
    gr_col='GRvalue'):
    pdata = median_std_data(separated_data,drug_col=drug_col, time_col=time_col, cell_col=cell_col, dose_col=dose_col, gr_col=gr_col)
    pdata = pdata.reset_index()
    pdata['condition'] = pdata[[drug_col, cell_col, dose_col, time_col]].apply(lambda x: '_'.join(x.astype(str)),axis = 1)
    return pdata

def waterfall_plot(data
    drug_col='agent', 
    time_col='time', 
    cell_col='cell_line', 
    dose_col='concentration', 
    gr_col='GRvalue',
    figsize=(16,9),
    color_sep='cell_line'):
    data.sort_values(dose_col,ascending=False,inplace=True)
    fig_width,fig_height = (16,9)
    plt.figure(figsize=figsize)
    width = fig_width/data.shape[0]
    num_colors = data[color_sep].unique().shape[0]
    for i, color in enumerate(data[color_sep].unique()):
        plot_data = data[data[color_sep]==color]
        sorted_dose_vector = sorted(data[dose_col].unique())
        c = ((i+1)/num_colors,0.5,1)
        color_table = pd.Series(sns.light_palette(c,input='hls',n_colors=len(sorted_dose_vector)),index = sorted_dose_vector)
        plt.bar(plot_data.index,plot_data.median_GR,yerr=plot_data.std_GR,width=width*2,color=color_table[plot_data.concentration],label=color)
    plt.title('DrugB', size=16)
    plt.legend(fontsize=16)
    plt.ylabel('GRvalue',size=16)
    plt.xlabel('Treatment conditions (dose shown as color saturation, low dose to low saturation)',size=16)
    plt.xticks([])