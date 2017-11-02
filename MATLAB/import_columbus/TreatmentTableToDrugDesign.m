function Design = TreatmentTableToDrugDesign(t_design, plate_dims)
% Design = TreatmentTableToDrugDesign(t_design, plate_size)
%   Convert a design table in an Design structure with standard fields for
%   treatment. Assume only one plate therefore enforces a single entry per
%   well.
%
%   t_design :  table with the following columns:
%                   - DrugName (and HMSLid)
%                   - Well
%                   - Conc (for concentration in uM)
%               and annotation columns:
%                   - DMSO (added to the treated_wells)
%               and optional columns (stored as 'Perturbations')
%                   - pert_type (e.g. trt_cp, vehicle_ctl, ...)
%                   - SeedingNumber
%                   - other perturbations (e.g. EGF/...)
%                   - DrugName2/3/... and Conc2/3/... for multiple drugs per well
%   plate_dims : specify the plate_size if not all well are treated
%
%   Design :    array of design structures with the following fields:
%                   - plate_dims (plate dimension)
%                   - treated_wells (wells treated with DMSO)
%                   - well_volume (in uL)
%                   - Drugs (structure with DrugName
%                       and layout - concentration given in uM)
%                   - Perturbations (structure with Name
%                       and layout - numeric array)
%


t_design = TableToCategorical(t_design,0);

assert(all(isvariable(t_design, {'DrugName' 'Well' 'Conc'})), ...
    'Need at least ''DrugName'' ''Well'' ''Conc'' as headers')

%% Check for additional drugs
Ndrugs = 1+max(isvariable(t_design, strcat('DrugName', cellfun(@(x) {num2str(x)}, num2cell(2:9)))));
assert(all(isvariable(t_design, strcat('Conc', cellfun(@(x) {num2str(x)}, num2cell(2:Ndrugs))))),...
    'Missing Conc2...%i for DrugName2...%i', Ndrugs, Ndrugs);

%% Check the wells and find the size of the plate
wells = AnyToString(t_design.Well);
assert(length(unique(wells)) == length(wells), ...
    'Multiple values for the same well found')

[rows, cols] = ConvertWellsToRowCol(wells);

if ~exist('plate_dims','var')
    plate_dims = 2.^ceil(log2([max(rows) max(cols)/1.5])).*[1 1.5];
end
assert((plate_dims(2)/plate_dims(1))==1.5)


%% Treated wells (all wells listed in the file)

treated_wells = false(plate_dims);
treated_wells(sub2ind(plate_dims, rows, cols)) = true;

if isvariable(t_design, 'DMSO')
    idx = t_design.DMSO>0;
    treated_wells(sub2ind(plate_dims, rows(idx), cols(idx))) = ...
        true | treated_wells(sub2ind(plate_dims, rows(idx), cols(idx)));
end

%% Drugs

DrugNames = setdiff(unique(t_design.DrugName), {'' '-'});
for iND2 = 2:Ndrugs
    DrugNames = unique([DrugNames; setdiff(unique(t_design.(['DrugName' num2str(iND2)])),...
        {'' '-'})], 'stable');
end
DrugNames = AnyToString(DrugNames);

% get the layout for each drug
layouts = cell(length(DrugNames),1);
for iD = 1:length(DrugNames)
    layouts{iD} = zeros(plate_dims);
    idx = t_design.DrugName==DrugNames{iD};
    layouts{iD}(sub2ind(plate_dims, rows(idx), cols(idx))) = t_design.Conc(idx);

    % case of multiple columns
    for iND2 = 2:Ndrugs
        idx = t_design.(['DrugName' num2str(iND2)])==DrugNames{iD};
        if any(layouts{iD}(sub2ind(plate_dims, rows(idx), cols(idx)))>0)
            warnprintf('%-14s found more than once in the same well --> summing', ...
                DrugNames{iD});
        end
        layouts{iD}(sub2ind(plate_dims, rows(idx), cols(idx))) = ...
            layouts{iD}(sub2ind(plate_dims, rows(idx), cols(idx))) + ...
            t_design.(['Conc' num2str(iND2)])(idx);
    end
end

% check for an HMSLid
if isvariable(t_design, 'HMSLid')
    HMSLids = cell(length(DrugNames),1);
    for iD = 1:length(DrugNames)
        HMSLids{iD} = AnyToString(t_design.HMSLid(find(t_design.DrugName==DrugNames{iD},1,'first')));
        if isempty(HMSLids{iD})
            % check for additional information in other columns
            for iND2 = 2:Ndrugs
                if isvariable(t_design, ['HMSLid' num2str(iND2)])
                HMSLids{iD} = AnyToString(t_design.(['HMSLid' ...
                    num2str(iND2)])(find(t_design.(['DrugName' num2str(iND2)])==...
                    DrugNames{iD},1,'first')));
                end
            end
        end
        if isempty(HMSLids{iD})
            HMSLids{iD} = '-';
        end
    end

    Drugs = struct('DrugName', DrugNames, 'layout', layouts, 'HMSLid', HMSLids);
else
    Drugs = struct('DrugName', DrugNames, 'layout', layouts);
end

%% Other perturbations

PertNames = ToColumn(setdiff(varnames(t_design), [{'Well' 'DrugName' 'Conc' 'Row' 'Cols' 'HMSLid'} ...
    strcat('Conc', cellfun(@(x) {num2str(x)}, num2cell(2:Ndrugs))) ...
    strcat('DrugName', cellfun(@(x) {num2str(x)}, num2cell(2:Ndrugs))) ...
    strcat('HMSLid', cellfun(@(x) {num2str(x)}, num2cell(2:Ndrugs)))], 'stable'));

layouts = cell(length(PertNames),1);
for iP = 1:length(PertNames)
    temp = t_design.(PertNames{iP});
    if iscategorical(temp);
        temp = cellstr(temp);
    end
    if isnumeric(temp)
        layouts{iP} = zeros(plate_dims);
    else
        layouts{iP} = repmat({''},plate_dims);
    end
    layouts{iP}(sub2ind(plate_dims, rows, cols)) = temp;
end
Perturbations = struct('Name',PertNames, 'layout', layouts);

fprintf(['\tAdded perturbations(s): ''' cellstr2str(PertNames', ''', ''') '''\n'])
%%

Design =  struct('plate_dims', plate_dims, 'treated_wells', treated_wells, ...
    'Perturbations', Perturbations, 'Drugs', Drugs);
