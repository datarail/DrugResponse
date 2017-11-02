%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% example of the D300/Columbus suite
%%%     - working from a hpdd file and modifying the layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the hpdd file and convert in a Design format
hpdd_filename = './results_20140620/20140620_SeedingAssay_4_MCF10A mCherry_HRG.hpdd';
[Design, temp_barcode] = hpdd_importer(hpdd_filename);

% add the perturbation to the Design
seeding = struct('Name','SeedingNumber', 'layout', ...
    (5e3.*(2.^-[5 5 0 0 1 1 2 2 3 3 4 4 5 5 0 0]'))*ones(1,24));
for i = 1:length(Design)
    Design(i).Perturbations(1) = seeding;
end

% save it
save ./results_20140620/20140620_SeedingAssay_4.mat Design

%% load and correct the barcode

t_plates = tsv2table('./results_20140620/20140620_barcode.tsv');

% adjuste some variable names
t_plates.Properties.VariableNames{'Replicate'} = 'DesignNumber';
t_plates.Properties.VariableNames{'Treatmentfile'} = 'TreatmentFile';
assert(all(eqtable(temp_barcode(:,{'TreatmentFile' 'DesignNumber'}), ...
    t_plates(:,{'TreatmentFile' 'DesignNumber'}))))
t_plates.Properties.VariableNames{'time_day_'} = 'Time';
t_plates.Time(:) = 72;

% change the treatment file to the new one
t_plates.TreatmentFile(:) = {'20140620_SeedingAssay_4.mat'};


%% load the data from Columbus

t_data = Import_PlatesCellCountData( ...
    './results_20140620/20140620_density effects_MCF10A mCherry_4 drugs with HRG.txt', ...
    t_plates);


%% Get only the data at the last time point (this was a time course)

t_data.Date = cellfun(@datenum, t_data.Date);
t_maxTP = collapse(t_data,@max, 'keyvars',{'Barcode'}, 'valvars', {'Date'});
t_data = t_data(ismember(t_data.Date, t_maxTP.Date),:);


%% Add the treatment to the data table

t_annotated = Annotate_CellCountData(t_data, './results_20140620/');
% correct the control, because the second drug is a marker
t_annotated.pert_type(t_annotated.Conc2==0) = {'ctl_vehicle'};


%% process the data (evaluate the relative cell count, do the mean across replicates)

[t_mean, t_processed] = Merge_CellCountData(t_annotated, {'SeedingNumber' 'modifier'});


%% fit the dose response curves and get the fit parameters

t_fits = ExtractCurves_CellCountData(t_mean, {'SeedingNumber' 'modifier'});


%% save the data

save DataFits_20140808.mat t_fits t_mean t_processed


%% plot some dose response curves

DrugNames = unique(t_processed.DrugName);

for i=1:length(DrugNames)
    get_newfigure(90+i,[10 10 1600 400], ['HRG_SeedingNumber_' char(DrugNames(i)) '.pdf'])

    plot_multidims(t_processed(t_processed.DrugName==DrugNames(i),:), 'yplotkey', 'CellLine', 'xplotkey', 'SeedingNumber', 'yspacing', .1, 'xaxiskey', ...
        'Conc', 'yaxiskey', 'RelCellCnt', 'xtransform', @log10, 'colorkey','modifier', ...
        'axischanges', @(x) set(x,'fontsize',6,'ylim',[0 1.5]))
end
