
function t_Tmean = TimeMerge_CellCountData(t_processed, NTimePlates, plate_inkeys, cond_inkeys, numericfields)
% t_Tmean = TimeMerge_CellCountData(t_processed, NTimePlates, plate_inkeys, cond_inkeys, numericfields)

%% assign and control the variables
plate_keys = {'CellLine' 'DeltaT'};
if exist('plate_inkeys','var') && ~isempty(plate_inkeys)
    plate_keys = unique([plate_keys plate_inkeys]);
end
plate_keys = intersect(plate_keys, varnames(t_processed));

cond_keys = [{'DrugName' 'Conc' 'Time' 'HMSLid' 'pert_type'} ...
    strcat('DrugName', num2cellstr(1:4)) strcat('Conc', num2cellstr(1:4))] ;
if exist('cond_inkeys','var') && ~isempty(cond_inkeys)
    cond_keys = unique([cond_keys cond_inkeys]);
end
cond_keys = intersect(cond_keys, varnames(t_processed));

Relvars = intersect({'RelCellCnt' 'RelGrowth' 'nRelGrowth'}, varnames(t_processed));
Relvars = Relvars(all(cellfun(@(x) isnumeric(x) & isscalar(x), table2cell(t_processed(1:min(40,end),Relvars)))));
labelfields = {'DesignNumber' 'Barcode' ...
    'Date' 'Row' 'Column' 'Well' 'TreatmentFile' 'Replicate'}; % remove all technical replicate info
if ~exist('numericfields','var')
    numericfields = setdiff(t_processed.Properties.VariableNames( ...
        all(cellfun(@(x) isnumeric(x) & isscalar(x), table2cell(t_processed(1:min(40,end),:))))), ...
        [plate_keys cond_keys labelfields Relvars]);
    if ~isempty(numericfields)
        fprintf('\tThese numeric fields will be averaged (set as cond_inkeys to use them as key):\n');
        for i=1:length(numericfields)
            fprintf('\t - %s\n', numericfields{i});
        end
    end
end

%%

% find the number of different of plates to merge and group them based on
% the plate_keys (with Time==0)
if isvariable(t_processed,'Time')
    t_processed.Time = round(t_processed.Time,2);
end
if isvariable(t_processed,'DeltaT')
    t_processed.DeltaT = round(t_processed.DeltaT,1);
end

t_plates = unique(t_processed(:,plate_keys));
t_Tmean = table;
var50 = varnames(t_processed);
var50 = var50(regexpcell(var50,'50')>0);
for iV = 1:length(var50)
    fprintf('  -> %s mean is performed in the log10 domain\n', var50{iV});
    t_processed.(var50{iV}) = log10(t_processed.(var50{iV}));
end


for iP = 1:height(t_plates)

    subt = t_processed(eqtable(t_processed, t_plates(iP,:)),:);

    Times = unique(subt.Time);
    Time = NaN*subt.Time;
    % find the best plate grouping
    mindiff = NaN(NTimePlates,1);
    for i=1:NTimePlates
        mindiff(i) = sum(sum(diff(reshape(Times(i:(i+...
            NTimePlates*floor((length(Times)-i)/NTimePlates)-1)),NTimePlates,[]))));
    end

    for iT = 0:ceil(length(Times)/NTimePlates)
        t = Times(max(1,argmin(mindiff)+NTimePlates*(iT-1)):...
            min(end,argmin(mindiff)+NTimePlates*iT-1));
        Time(ismember(subt.Time, t)) = mean(t);
    end
    Time = round(Time,2);
    temp = [table(Time) subt(:,setdiff(varnames(subt),'Time','stable'))];

    otherkeys = setdiff(varnames(temp), [cond_keys plate_keys], 'stable');
    otherkeys = otherkeys(all(cellfun(@(x) (isnumeric(x) | ischar(x) | iscategorical(x)) & length(x)==1, ...
        table2cell(temp(1:min(10,end),otherkeys)))));
    otherkeys = otherkeys(colfun_array(@(x) height(unique(x))==1, ...
        temp(:,otherkeys))==1);
    t_Tmean = [t_Tmean; collapse(temp, @mean, 'keyvars', setdiff([cond_keys plate_keys otherkeys], ...
        [labelfields Relvars numericfields],'stable'), 'valvars', [Relvars numericfields])];
end


for iV = 1:length(var50)
    t_Tmean.(var50{iV}) = 10.^(t_Tmean.(var50{iV}));
end

% round the DeltaT to avoid artifacts
if isvariable(t_Tmean, 'DeltaT')
    t_Tmean.DeltaT = round(t_Tmean.DeltaT,2);
    t_Tmean.T0 = round(t_Tmean.Time - .5*t_Tmean.DeltaT,2);
    t_Tmean.Tend = round(t_Tmean.Time + .5*t_Tmean.DeltaT,2);
end
t_Tmean = sortrows(t_Tmean, plate_keys);
