function Designs = TreatmentDesign(DrugNames, HMSLids, SingleDoses, nReps, varargin)
% Designs = TreatmentDesign(DrugNames, HMSLids, SingleDoses, nReps, varargin)
%
%   Produce randomized treatment desgins based on the specifications
%
%   DrugNames       label used to identify the drug
%   HMSLid          corresponding HMSL id
%   SingleDoses     cell array of doses vectors, order corresponds to
%                       DrugNames, dispensed as a single agent
%   nReps           Number of repeats (randomized replicates)
%
%
%   Optional arguments (name-value pairs):
%
%   ComboDoses      cell array of doses vectors, order corresponds to
%                       DrugNames
%   DrugPairs       [Nx2] list of drug combinations to test (with
%                       ComboDoses)
%   Seed            initial seed, set to 0 for no randomization (default 1).
%                       For each repeat the seed is 'Seed+i-1'
%   edge_ctrl       default is true
%   stock_conc      in uM, default is 1e4 (= 10mM)
%   Vehicle         drug vehicle (default is DMSO)
%   well_volume     in uL, default is 60 (for 384-well plates)
%   plate_dims      either [N of rows ,  N of columns] or [N of wells],
%                       default is [16 24]
%   max_DMSOpct     maximum percent of DMSO dispensed (default 0.2%)
%   min_volume      step volume dispensed (in uL, default is 2e-5)
%   step_volume     step volume dispensed (in uL, default is 1.3e-5)
%   Perturbations   add the preturbations directly in the structure (same
%                       for all replicates)
%
%   Output: array od Design structure matching the DrugResponseSuite
%
%


%% inputs and parser
p = inputParser;

addParameter(p,'ComboDoses',[], @(x) all(cellfun(@isnumeric,x)));
addParameter(p,'ComboLists',[], @(x) all(cellfun(@isnumeric,x(:))));
addParameter(p,'DrugPairs',zeros(0,2), @isnumeric);
addParameter(p,'Seed',1,  @(x) isscalar(x) && isnumeric(x));
addParameter(p,'edge_ctrl',true, @(x) islogical(x) & isscalar(x));
addParameter(p,'stock_conc',1e4, @(x) isvector(x) && isnumeric(x));     % in uM
addParameter(p,'Vehicle','DMSO', @(x) iscellstr(x) || ischar(x));     % string
addParameter(p,'well_volume',60, @(x) isscalar(x) && isnumeric(x));       % in uL
addParameter(p,'plate_dims',[16 24], @(x) isvector(x) && isnumeric(x));
addParameter(p,'Perturbations',repmat(struct('Name',[],'layout',[]),0,0), @isstruct);

% based on the specifications of the D300
addParameter(p,'min_volume', 1.3e-5, @(x) isscalar(x) && isnumeric(x)); % minimum volume of 13pl (in uL)
addParameter(p,'step_volume', 2e-5, @(x) isscalar(x) && isnumeric(x)); % minimum step of 20pl (in uL)
addParameter(p,'max_DMSOpct', .2, @(x) isscalar(x) && isnumeric(x)); % maximum 0.2% percent of DMSO

parse(p, varargin{:})
p = p.Results;

% avoid too high level of DMSO (max is 0.2%)
max_volume = p.well_volume*p.max_DMSOpct/100;

%% checks
DrugNames = ToColumn(DrugNames);
assert((length(DrugNames)==length(HMSLids)) || isempty(HMSLids))
assert(length(DrugNames)==length(SingleDoses))

assert(isempty(p.DrugPairs) || (~isempty(p.ComboDoses) || ~isempty(p.ComboLists)))
assert(all(p.DrugPairs(:)<=length(DrugNames)) || isempty(p.DrugPairs))
assert((length(DrugNames)==length(p.ComboDoses)) || isempty(p.ComboDoses))
assert(size(p.DrugPairs,2)==2)
assert(all(size(p.DrugPairs)==size(p.ComboLists)) || isempty(p.ComboLists))

stock_conc = p.stock_conc;
assert((length(DrugNames)==length(stock_conc)) || (length(stock_conc)==1))
assert(all(stock_conc>20 & stock_conc<5e4), 'Stock concentration should be in uM')
if length(stock_conc)==1
    stock_conc = stock_conc*ones(length(DrugNames),1);
end
stock_conc = num2cell(ToColumn(stock_conc));

Vehicle = p.Vehicle;
assert((iscellstr(Vehicle) && length(DrugNames)==length(Vehicle)) || ...
    (ischar(Vehicle)))
if ischar(Vehicle)
    Vehicle = repmat({Vehicle},length(DrugNames),1);
end

if length(p.plate_dims)==1
    p.plate_dims = sqrt(p.plate_dims*[1/1.5 1.5]);
end
warnassert(round(log2(p.plate_dims(1)))==log2(p.plate_dims(1)), ...
    'Non-standard plate size; edge_ctrl may be not working properly')

nWells = prod(p.plate_dims);
if p.well_volume<1e2 && nWells<100
    warnprintf('Default volume is %.0f set for 384-well plate; check if correct', ...
        p.well_volume)
end
assert(p.well_volume>1 && p.well_volume<1e4, 'Well volume should be given in uL')

if p.Seed==0 && nReps>1
    warnprintf('Seed is set to 0 (i.e. no randomization) and Nreps>1')
    warnprintf('This means all replicates will have the same layout!')
    pause(5)
end

for iP=1:length(p.Perturbations)
    if ~isempty(p.Perturbations(iP).Name)
        assert(all(size(p.Perturbations(iP).layout)==p.plate_dims), ...
            'Plate dims do not match Perturbation (%s)', p.Perturbations(iP).Name)
    end
end

%% initialization


Drugs = struct('DrugName', DrugNames, 'HMSLid', ToColumn(HMSLids),...
    'stock_conc', stock_conc, 'layout', zeros(p.plate_dims), 'Vehicle', ToColumn(Vehicle));
Designs = struct('plate_dims', repmat({p.plate_dims}, nReps, 1), ...
    'treated_wells', repmat({true(p.plate_dims)}, nReps, 1), ...
    'Vehicle', repmat({repmat({''},p.plate_dims)}, nReps, 1), ...
    'well_volume', repmat({p.well_volume}, nReps, 1), ...
    'Drugs', repmat({Drugs}, nReps, 1), 'Seed', num2cell(p.Seed-1+(1:nReps)'), ...
    'Perturbations', repmat({p.Perturbations}, nReps, 1));

SingleDoses = cellfun2(@(x,y,z) round_Doses(x,y,z,'Single',p.min_volume, ...
    p.step_volume, max_volume, p.well_volume), SingleDoses, stock_conc, DrugNames);

if isempty(p.ComboLists) && ~isempty(p.ComboDoses)
    p.ComboLists = [p.ComboDoses(p.DrugPairs(:,1)) p.ComboDoses(p.DrugPairs(:,2))];
end

if ~isempty(p.ComboLists)
    p.ComboLists(:) = cellfun2(@(x,y,z) round_Doses(x,y,z,'Combo',p.min_volume, ...
        p.step_volume, max_volume, p.well_volume), p.ComboLists(:), ...
        stock_conc(p.DrugPairs(:)), DrugNames(p.DrugPairs(:)));
    fprintf('\n')
    for iD = 1:length(SingleDoses)
        if any(~ismember([p.ComboLists{p.DrugPairs==iD}], SingleDoses{iD}))
            warnprintf('Some doses for the combo are not part of the single doses for %s', ...
                DrugNames{iD})
        end
    end
    nTreatments = sum(cellfun(@length,SingleDoses)) + ...
        sum(cellfun(@(x,y) length(x)*length(y), p.ComboLists(:,1), ...
        p.ComboLists(:,2)));
else
    nTreatments = sum(cellfun(@length,SingleDoses));
end



%% define the position of the fixed controls
fprintf('Number of test wells: %i\n',nTreatments);
ctrl_cnt = nWells - nTreatments;
if ctrl_cnt<0
    error('Too many well used (%i out of %i)', nTreatments, nWells)
elseif ctrl_cnt<6
    warnprintf('Too many well used (%i out of %i), need at least 6 control wells',...
        nTreatments, nWells)
    fprintf('Confirm:')
    pause
    fprintf('\n')
end



if p.edge_ctrl
    [ctrlidx, treated_wells] = DefineFixedControlPositions(p.plate_dims, ctrl_cnt);
else
    ctrlidx = [];
    treated_wells = true(p.plate_dims);
end
for iR=1:nReps
    Designs(iR).treated_wells = treated_wells;
end
%% define all possible treatments
allTreatments = zeros(length(DrugNames), nTreatments);
cnt = 0;
for iD = 1:length(DrugNames)
    allTreatments(iD, cnt+(1:length(SingleDoses{iD}))) = SingleDoses{iD};
    cnt = cnt + length(SingleDoses{iD});
end
for iCo = 1:size(p.DrugPairs,1)
    for iD1 = 1:length(p.ComboLists{iCo,1})
        allTreatments(p.DrugPairs(iCo,1), cnt+(1:length(p.ComboLists{iCo,2}))) = ...
            p.ComboLists{iCo,1}(iD1);
        allTreatments(p.DrugPairs(iCo,2), cnt+(1:length(p.ComboLists{iCo,2}))) = ...
            p.ComboLists{iCo,2};
        cnt = cnt + length(p.ComboLists{iCo,2});
    end
end

assert(cnt==nTreatments)
assert(all(any(allTreatments>0)))


%% assign the design for each replicate
for iR = 1:nReps
    
    Designs(iR).Drugs = RandomizePlatePositions(Designs(iR).Drugs, ...
        allTreatments, nWells, p.plate_dims, ctrlidx, ctrl_cnt, (p.Seed~=0)*(p.Seed+iR-1));
    
    allDrugs = reshape([Designs(iR).Drugs.layout], [p.plate_dims length(DrugNames)]);
    nDrugs = sum(allDrugs>0,3);
    assert(all(squeeze(sum(sum((allDrugs>0).*repmat(nDrugs==1,1,1,length(DrugNames)),2),1))==...
        cellfun(@length,SingleDoses)))
    for iCo = 1:size(p.DrugPairs,1)
        assert(sum(sum( all(allDrugs(:,:,p.DrugPairs(iCo,:))>0,3)))== ...
            (length(p.ComboLists{iCo,1})*length(p.ComboLists{iCo,2})*...
            sum(all(p.DrugPairs==(ones(size(p.DrugPairs,1),1)*p.DrugPairs(iCo,:)),2))))
    end
    
    Designs(iR) = SetDrugVehicle(Designs(iR));
    
end


end

