function t_rate = TimeCourse_DivRate(t_data, cond_inkeys)
% t_rate = TimeCourse_DivRate(t_data, cond_inkeys)
%
%   process the data from cell count. The data will be split for controls according
%   to plate_keys (with 'Barcode' 'CellLine' 'Time' as mandatory and
%   default).
%   Needs also the column 'pert_type' (looking for value
%   'ctl_vehicle' or 'trt_cp'/'trt_poscon'), 'Barcode' and 'Well' as unique
%   identifiers
%
%   outputs are:    t_rate with the values for all treated wells
%
%   NOTE: not handling a 'real' day 0 --> need to be improved!

% unique identifiers for following timecourse
if ~exist('cond_inkeys', 'var')
    trace_vars = {};
else
    trace_vars = cond_inkeys;
end
for iF = ['Barcode' 'Well' 'DrugName' 'Conc' 'CellLine' ...
        strcat('DrugName',num2cell('2':'9')) strcat('Conc',num2cell('2':'9'))]
    if isvariable(t_data, iF{:})
        trace_vars = [trace_vars iF];
    end
end
trace_vars = unique(trace_vars);

if exist('plate_inkeys','var') && ~isempty(cond_inkeys)
    plate_keys = unique([{'CellLine' 'Barcode' 'Time'} cond_inkeys]);
else
    plate_keys = {'CellLine' 'Barcode' 'Time'};
end
plate_keys = intersect(plate_keys, varnames(t_data), 'stable');
%%

t_location = unique(t_data(:,trace_vars));

t_rate = table();
for it = 1:height(t_location)
    % get all time points for each well
    t_temp = sortrows(t_data(eqtable(t_data(:,trace_vars), t_location(it,:)),:),'Time');

    t_temp = [t_temp(1,:); t_temp];
    t_temp.Time(1) = 0;
    if isvariable(t_temp, 'Day0Cnt')
        t_temp.Cellcount(1) = t_temp.Day0Cnt(1);
    else
        t_temp.Cellcount(1) = t_temp.Cellcount(2);
    end
    
    Time = mean([t_temp.Time(1:(end-1)) t_temp.Time(2:end)],2);
    Cellcount =  mean([t_temp.Cellcount(1:(end-1)) t_temp.Cellcount(2:end)],2);
    dx = diff(t_temp.Cellcount);
    dt = diff(t_temp.Time)/24;
    DivRate = dx./dt./Cellcount;
    DivRate = smooth(DivRate, 3);
    
    n = NaN(width(t_temp),1);
    for i=1:width(t_temp), n(i) = length(unique(t_temp.(i))); end;
    annotation_vars = setdiff(t_temp.Properties.VariableNames(n==1), ...
        [trace_vars 'RelCellCnt' 'RelGrowth' 'Day0Cnt' 'Cellcount']);
    if ~isempty(t_rate)
        annotation_vars = intersect(varnames(t_rate), annotation_vars,'stable');
    end
    t_rate = [t_rate;
        t_temp(2:end, trace_vars) table(Time, Cellcount, DivRate) ...
        t_temp(2:end,annotation_vars)];
end

%
if isvariable(t_rate, 'pert_type')
    t_ctrl = collapse(t_rate(t_rate.pert_type=='ctl_vehicle', ['DivRate' plate_keys]), ...
        @(x)mean(max(x,0)), 'keyvars', plate_keys);
    
    t_rate.RelDivRate = NaN(height(t_rate),1);
    for i = 1:height(t_ctrl)
        idx = eqtable(t_ctrl(i,plate_keys), t_rate(:,plate_keys));
        t_rate.RelDivRate( idx ) = t_rate.DivRate( idx )./t_ctrl.DivRate(i);
    end
    t_rate.RelDivRate(~isnan(t_rate.RelDivRate)) = ...
        max(min(t_rate.RelDivRate(~isnan(t_rate.RelDivRate)), 4), -2);
end
