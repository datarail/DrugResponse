dates = {'20140625' '20140722' '20140826'};
filenames = {'./Cellcount/20140625_LJP_six cell lines_Hoechst.txt', ...
    './Cellcount/20140722_LJP_six cell lines_Hoechst LDR.txt', ...
    './Cellcount/20140826_LJP_six cell lines_Hoechst.txt'};

for id = 1:length(dates)
    %%
    t_data = Import_PlatesCellCountData(filenames{id}, ...
        ['barcode_' dates{id} '.tsv']);
    %%

    t_annotated = Annotate_CellCountData(t_data);

    %%

    [t_mean, t_processed] = Merge_CellCountData(t_annotated, [], {'plate'}););

    %%

    save(['data_' dates{id} '.mat'], 't_annotated', 't_mean', 't_processed')
end
