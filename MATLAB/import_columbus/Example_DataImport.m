

fprintf('\n\nExample of importing the data:\n---------------------------\n\n');

fprintf('folder = ''./results/''; \n')

fprintf('t_data = Import_PlatesCellCountData([folder ''ColumbusOutput.txt''], ... \n');
fprintf('\t[folder ''PlateInfo.tsv''], ''NobjField'', ... %% example with multiple non-default fields \n');
fprintf('\t{''Nuclei_Hoechst_NumberOfObjects'' ''Nuclei_LDR_NumberOfObjects'' ...\n')
fprintf('\t\t''Nuclei_Hoechst_LDRpos_NumberOfObjects''}, ''Cellcount'', @(x) x(:,1)-x(:,3) );\n');

fprintf('\n%%%%\nt_annotated = Annotate_CellCountData(t_data, folder); %% where the treatment design files are stored \n');
fprintf('t_corrected = EdgeCorrecting_CellCountData(t_annotated); %% if edge correction is needed \n');
fprintf('TestBias_multiplates(t_corrected); %% test if there is a bias\n');

%
fprintf('\n%%%%\n[t_mean, t_processed] = Merge_CellCountData(t_corrected); \n');

fprintf('\n%%%%\nt_fits = ExtractCurves_CellCountData(t_mean, keys, pcutoff); \n');
