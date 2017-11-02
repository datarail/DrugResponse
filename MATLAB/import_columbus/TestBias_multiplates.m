function [BiasValue, BiasResults] = TestBias_multiplates(t_data, BiasCutoff, plotting, valvars)
% [BiasValue, BiasResults] = TestBias_multiplates(t_data, BiasCutoff, plotting, valvars)
% t_data:     need variables 'Well' or 'Row'/'Column'; all other
%                         columns will be checked against.
%
% BiasCutoff: [minimal drop cutoff      p-value cutoff];
%                         default = [.1 .01]
%                     can be a 3x2 matrix for each condition:
%                             edge, column, row
%
%           plotting :  = 1 , pause for every plot
%                       = .5, pause 1s for the biased ones
%                       = 2 (or string for file name), save everything as a pdf
%


if ~exist('BiasCutoff','var') || isempty(BiasCutoff)
    BiasCutoff = [1;1;1]*[.1 .01];
elseif all(size(BiasCutoff)==[1 2])
    BiasCutoff = [1;1;1]*BiasCutoff;
elseif ~all(size(BiasCutoff)==[3 2])
    error('Wrong size for BiasCutoff input, needs to be 1x2 or 3x2')
end


if ~exist('plotting','var') || isempty(plotting)
    plotting = .5;
end
if ischar(plotting)
    FileName = [plotting '.pdf'];
    plotting = 2;
elseif plotting == 2
    FileName = 'TestBias_results.pdf';
end

if ~exist('valvars','var') || isempty(valvars)
    valvars = 'Cellcount';
end

t_plates = unique(t_data(:,{'Barcode' 'Time'}));

BiasValue = zeros(height(t_plates), 3);
clear BiasResults;

get_newfigure(999, [40 100 600 400]);
if plotting==2
    mkdir temp_pdf
end

for ip = 1:height(t_plates)
    t_plate = t_data(eqtable(t_data,t_plates(ip,:)),intersect(varnames(t_data), ...
        [{'Barcode' 'Time' 'Column' 'Row' 'Well'} valvars]));

    bias_res = cell(1,3);
    [biased, bias_res{:}] = TestPlateBias(t_plate, BiasCutoff, plotting);

    BiasResults(ip) = struct('PlateInfo', t_plates(ip,:), 'edge_res', bias_res{1},...
        'col_res', bias_res{2}, 'row_res', bias_res{3});

    if plotting == 2
        savegcf(['./temp_pdf/fig_' num2str(ip,'%03i') '.pdf'])
    end

    for i=find(biased)
        BiasValue(ip,i) = min(BiasValue(ip,i), min(bias_res{i}(...
            bias_res{i}(:,1)<-BiasCutoff(i,1) & bias_res{i}(:,3)<BiasCutoff(i,2),1)));
    end
end

if plotting == 2
    merge_pdf('./temp_pdf/*pdf', FileName , true)
    pause(.3)
    rmdir temp_pdf
end
