function t_fits = ExtractCurves_CellCountData(t_data, keys, pcutoff)
% t_fits = ExtractCurves_CellCountData(t_data, keys, pcutoff)
%   Sigmoidal fit on the drug response data (expect concentration in uM) to
%   extract the following parameters:
%   If the relative cell count (column RelCellCnt) is measured:
%       - IC50
%       - Einf (effect -> = 0 if no drug response)
%       - Hill
%       - Area
%       - EC50
%   If the relative growth (column RelGrowth) is measured:
%       - GI50
%       - GIinf
%   If the normalized relative growth (column GRvalue) is measured:
%       - GR50
%       - GRinf
%       - Hill
%       - Area
%       - EC50
%   All the outputs are saved in a table with annotations including drug
%   concentrations and initial values.
%
%   keys are used for aggregation of the data; default are : CellLine,
%   DrugName, Time, SeedingNumber, Date.
%
%   pcutoff is the cutoff for the p-value of a F-test against a flat line.
%


%%
if exist('keys','var')
    keys = intersect(t_data.Properties.VariableNames, ...
        [{'CellLine' 'DrugName' 'Time' 'Date' 'SeedingNumber' 'SeedingDensity'} keys]);
else
    keys = intersect(t_data.Properties.VariableNames, ...
        {'CellLine' 'DrugName' 'Time' 'SeedingNumber' 'Date' 'SeedingDensity'});
end

if ~isempty(setdiff(intersect(varnames(t_data), ...
        strcat('DrugName', cellfun(@(x) {num2str(x)}, num2cell(2:9)))), ...
        keys))
    warnprintf('Mutiple drugs found in the data, but only the first one is a key; removing any double agent')
    MultiDrugs = setdiff(intersect(varnames(t_data), ...
        strcat('DrugName', cellfun(@(x) {num2str(x)}, num2cell(2:9)))), ...
        keys);
    for i=1:length(MultiDrugs)
        t_data(t_data.(MultiDrugs{i})~='-',:) = [];
    end
end

if ~exist('pcutoff','var')
    pcutoff = .05;
end
fitopt.pcutoff = pcutoff;

fitopt2 = fitopt;

t_keys = unique(t_data(:,keys));


t_fits = table;
for ik = 1:height(t_keys)
    loop_waitbar(ik, height(t_keys))
    %%
    subt = sortrows(t_data(eqtable(t_keys(ik,:), t_data(:,keys)),:),'Conc');
    % temporary hack for compatibility
    if ismember('nRelGrowth', subt.Properties.VariableNames)
        subt.GRvalue = subt.nRelGrowth;
    end

    if height(subt)<4
        warnprintf(['Not enough data point for ' strjoin(table2cellstr(t_keys(ik,:),0))])
        continue
    end

    if all(subt.Conc<1e-3) || all(subt.Conc>1e3)
        warnprintf('Concentrations are expected in uM; fitopt have constraints')
    end

    t_temp = t_keys(ik,:);

    if ismember('RelCellCnt', subt.Properties.VariableNames)
        [IC50, Hill, Einf, Emax, AUC, r2, EC50, fit] = ...
            ICcurve_fit(subt.Conc, subt.RelCellCnt, 'IC50', fitopt);

        t_temp = [t_temp table(IC50, Hill, Einf, Emax, AUC, r2, EC50) ...
            table({fit}, {subt.Conc'}, {subt.RelCellCnt'}, 'VariableNames', ...
            {'fit' 'Conc' 'RelCellCnt'})];
    end

    if ismember('RelGrowth', subt.Properties.VariableNames)
        if all(ismember({'Day0Cnt' 'Ctrlcount'}, subt.Properties.VariableNames));
            fitopt2.ranges = [
                min(min(-subt.Day0Cnt./(subt.Ctrlcount-subt.Day0Cnt))*1.1,-.01) 1    %GImax
                max(min(subt.Conc)*1e-3,1e-7) min(max(subt.Conc)*1e2, 1e3)  %E50
                .1 5    % HS
                ]';
        else
            fitopt2 = fitopt;
        end
        [GI50, ~, GIinf, GImax, GIArea, GI_r2, ~, GI_fit] = ...
            ICcurve_fit(subt.Conc, subt.RelGrowth, 'GI50', fitopt2);
        t_temp = [t_temp table(GI50, GIinf, GImax, GIArea, GI_r2) ...
            table({GI_fit}, {subt.RelGrowth'}, 'VariableNames', {'GI_fit' 'RelGrowth'})];
    end

    if ismember('GRvalue', subt.Properties.VariableNames)
        [GR50, Hill_GR, GRinf, GRmax, GR_AOC, GR_r2, GEC50, GR_fit] = ...
            ICcurve_fit(subt.Conc, subt.GRvalue, 'GR50', fitopt);
        t_temp = [t_temp table(GR50, Hill_GR, GRinf, GRmax, GR_AOC, GEC50, GR_r2) ...
            table({GR_fit}, {subt.GRvalue'}, 'VariableNames', {'GR_fit' 'GRvalue'})];
    end

t_fits = [t_fits; t_temp];

end
