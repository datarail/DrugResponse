function t_fits = ExtractCurves_CellCountData2(t_data, keys, pcutoff)
% t_fits = ExtractCurves_CellCountData2(t_data, keys, pcutoff)
%   Sigmoidal fit on the drug response data (expect concentration in uM) to
%   extract the following parameters based on the GR_s and GR_d values
%       - GRs50
%       - GRsinf
%       - GRsmax
%       - h_GRs
%       - GRs_AOC
%       - GECs50
%       - GRd50
%       - GRdinf
%       - GRdmax
%       - h_GRd
%       - GRd_AOC
%       - GECd50
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
    
    if height(subt)<4
        warnprintf(['Not enough data point for ' strjoin(table2cellstr(t_keys(ik,:),0))])
        continue
    end

    if all(subt.Conc<1e-3) || all(subt.Conc>1e3)
        warnprintf('Concentrations are expected in uM; fitopt have constraints')
    end

    t_temp = t_keys(ik,:);

    if ismember('GR_s', subt.Properties.VariableNames)
        [GRs50, h_GRs, GRsinf, GRsmax, GRs_AOC, GRs_r2, GECs50, GRs_fit] = ...
            ICcurve_fit(subt.Conc, subt.GR_s, 'GR50', fitopt);
        t_temp = [t_temp table(GRs50, GRsmax, GRs_AOC, GRsinf, h_GRs, GECs50, GRs_r2) ...
            table({GRs_fit}, {subt.GR_s'}, 'VariableNames', {'GRs_fit' 'GR_s'})];
    end

    if ismember('GR_d', subt.Properties.VariableNames)
        [GRd50, h_GRd, GRdinf, GRdmax, GRd_AOC, GRd_r2, GECd50, GRd_fit] = ...
            ICcurve_fit(subt.Conc, subt.GR_d+1, 'GR50', fitopt);
        GRdinf = GRdinf-1; GRdmax = GRdmax-1; 
        GRd_fit = @(x) GRd_fit(x)-1;
        t_temp = [t_temp table(GRd50, GRdmax, GRd_AOC, GRdinf, h_GRd, GECd50, GRd_r2) ...
            table({GRd_fit}, {subt.GR_d'}, 'VariableNames', {'GRd_fit' 'GR_d'})];
    end
    
    if ismember('GRvalue', subt.Properties.VariableNames)
        [GR50, h_GR, GRinf, GRmax, GR_AOC, GR_r2, GEC50, GR_fit] = ...
            ICcurve_fit(subt.Conc, subt.GRvalue, 'GR50', fitopt);
        t_temp = [t_temp table(GR50, GRmax, GR_AOC, GRinf, h_GR, GEC50, GR_r2) ...
            table({GR_fit}, {subt.GRvalue'}, 'VariableNames', {'GR_fit' 'GR'})];
    end

t_fits = [t_fits; t_temp];

end
