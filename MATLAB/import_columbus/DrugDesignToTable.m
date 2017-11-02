function t_design = DrugDesignToTable(Design1, Perturbations, DrugOrder)
% t_design = DrugDesignToTable(Design, Perturbations, DrugOrder)
%   Convert a Design structure (not an array) with standard fields for
%   treatment into a design table. Order of Perturbations can be enforced
%   as well as Priority for the assignment of drugs in DrugName1/2/..
%
%
%   Design :    array of design structures with the following fields:
%                   - plate_dims (plate dimension)
%                   - treated_wells (wells treated with DMSO/drug/perturbation)
%                   - Drugs (structure with DrugName
%                       and layout - concentration given in uM)
%                   - Perturbations (structure with Name
%                       and layout - numeric array)
%
%   Perturbations:  list of Perturbations found in Design
%
%   DrugOrder : list of drugs to prioritize for storing in DrugName,
%                   DrugName2/... ; useful for comparison in designs with
%                   multiple drugs per well.
%
%   t_design :  table with the following columns:
%                   - Well
%                   - DrugName (and HMSLid)
%                   - Conc (for concentration in uM)
%               and annotation columns:
%                   - DMSO (added to the treated_wells)
%               and optional columns (stored as 'Perturbations')
%                   - pert_type (e.g. trt_cp, vehicle_ctl, ...)
%                   - SeedingNumber
%                   - other perturbations (e.g. EGF/...)
%                   - DrugName2/3/... and Conc2/3/... for multiple drugs per well
%
%               add the 'pert_type' if not present in Design1
%


%%
if isfield(Design1, 'Perturbations')
    trtwells = Design1.treated_wells ;
    untrtwells = ~Design1.treated_wells;
    for i=1:length(Design1.Perturbations)
        if iscell(Design1.Perturbations(i).layout)
            trtwells = trtwells | ~cellfun(@isempty,Design1.Perturbations(i).layout);
            untrtwells = untrtwells & ~cellfun(@isempty,Design1.Perturbations(i).layout);
        else
            trtwells = trtwells | Design1.Perturbations(i).layout>0;
            untrtwells = untrtwells & Design1.Perturbations(i).layout>0;
        end
    end
    [rows, cols] = find(trtwells);
    [Untrtrows, Untrtcols] = find(untrtwells);
    UntrtWell = ConvertRowColToWells(Untrtrows, Untrtcols);
else
    [rows, cols] = find(Design1.treated_wells);
    UntrtWell = {};
end
Well = ConvertRowColToWells(rows, cols);

% get all drugs and match the order
if ~isfield(Design1,'Drugs') || isempty(Design1.Drugs)
    DrugNames = {};
else
    DrugNames = {Design1.Drugs.DrugName};
    if exist('DrugOrder','var')
        [temp, order] = ismember(DrugOrder, DrugNames);
        order = [order(temp) find(~ismember(DrugNames, DrugOrder))];
        DrugNames = DrugNames(order);
    else
        order = 1:length(DrugNames);
    end
end
% get all perturbations
if isfield(Design1, 'Perturbations')
    PertNames = {Design1.Perturbations.Name};
    if exist('Perturbations','var') && ~isempty(Perturbations)
        PertNames = intersect(Perturbations, PertNames, 'stable');
        if isempty(PertNames)
            warnprintf('No perturbation found in Design that matches ''Perturbations''');
        end
    end
    fprintf('\tAdding perturbations: %s\n', strjoin(PertNames, ', '));
else
    PertNames = {};
end

if isfield(Design1.Drugs,'HMSLid') && ~isempty(Design1.Drugs)
    HMSLids = {Design1.Drugs.HMSLid};
    t_HMSLids = table([DrugNames {'-'}]', [HMSLids(order) {'-'}]', ...
        'VariableNames', {'DrugName' 'HMSLid'});
    t_HMSLids.HMSLid(cellfun(@isempty,t_HMSLids.HMSLid)) = {'-'};
else
    %%%% to be replaced by a link to the LINCS database for good
    %%%% annotations
    t_HMSLids = table([DrugNames {'-'}]', repmat({'-'},length(DrugNames)+1,1), ...
        'VariableNames', {'DrugName' 'HMSLid'});
end

%%

% determine the number of Drug Columns
if ~isfield(Design1,'Drugs') || isempty(Design1.Drugs)
    Ndrugs = 1;
    DrugConc = zeros([Design1.plate_dims 0]);
else
    Ndrugs = 1;
    DrugConc = reshape([Design1.Drugs(order).layout], ...
        [Design1.plate_dims length(DrugNames)]);
    if any(any(sum(DrugConc>0,3)>Ndrugs))
        Ndrugs = max(max(sum(DrugConc>0,3)));
        warnprintf('some wells have %i drugs, additional columns in output', ...
            Ndrugs)
    end
    % this is not the optimal way of storing multiple drugs because of the
    % hierarcy between DrugName and Conc as well as the redudancy and possible
    % swapping between Drug1 and Drug2 ; it makes matching between condition hard
end


% set the table columns
Conc = zeros(length(Well),1);
DrugName = repmat({'-'}, length(Well),1);
t_design = table(Well, DrugName, Conc);
for iAD = 2:Ndrugs
    t_design = [t_design table(DrugName, Conc, 'VariableNames', ...
        {sprintf('DrugName%i', iAD), sprintf('Conc%i', iAD)})];
end


% construct the drug columns
for iW = 1:height(t_design)
    Didx = find(DrugConc(rows(iW), cols(iW),:));
    if isempty(Didx), continue, end

    for iD = 1:length(Didx)
        t_design.(sprintf('DrugName%i', iD(iD>1))){iW} = DrugNames{Didx(iD)};
        t_design.(sprintf('Conc%i', iD(iD>1)))(iW) =  DrugConc(rows(iW), cols(iW),Didx(iD));
    end
end

% add the HMSLid
for iD = 1:Ndrugs
    temp = join(t_design(:,sprintf('DrugName%i', iD(iD>1))), t_HMSLids, ...
        'leftkeys', sprintf('DrugName%i', iD(iD>1)), ...
        'rightkeys', 'DrugName', 'rightvariables', 'HMSLid');
    idx = find(strcmp(varnames(t_design), sprintf('DrugName%i', iD(iD>1))));
    t_design = [t_design(:,1:idx) varnames(temp(:,'HMSLid'),{sprintf('HMSLid%i', iD)}), ...
        t_design(:,(idx+1):end)];
end
% correct the small trick to avoid having the same column name
t_design.Properties.VariableNames{'HMSLid1'} = 'HMSLid';

% add the perturbation columns
for iP = 1:length(PertNames)
    Pertvals = Design1.Perturbations(strcmp({Design1.Perturbations.Name}, ...
        PertNames{iP})).layout(sub2ind(Design1.plate_dims, rows, cols));

    t_design = [t_design, table(Pertvals, 'VariableNames', PertNames(iP))];
end

if ~isvariable(t_design, 'pert_type')
    pert_type = repmat({'trt_cp'}, height(t_design),1);
    pert_type(t_design.Conc==0 & ~ismember(t_design.Well,UntrtWell)) = {'ctl_vehicle'};
    pert_type(t_design.Conc==0 & ismember(t_design.Well,UntrtWell)) = {'Untrt'};
    t_design = [t_design table(pert_type)];
end

% add the vehicles
if isfield(Design1, 'Vehicle')
    Vehicles = Design1.Vehicle(sub2ind(Design1.plate_dims, rows, cols));
    Vehicles(cellfun(@isempty, Vehicles)) = {'-'};
    t_design = [t_design, table(Vehicles, 'VariableNames', {'Vehicle'})];
end

t_design = TableToCategorical(t_design, 0);
