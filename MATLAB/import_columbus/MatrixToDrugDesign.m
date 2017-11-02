function Designs = MatrixToDrugDesign(DrugNames, HMSLids, DrugLayout, nReps, varargin)
% Designs = MatrixToDrugDesign(DrugNames, HMSLids, DrugLayout, nReps, varargin)
%
%   Produce randomized treatment desgins based on the specifications
%
%   DrugNames       label used to identify the drug
%   HMSLid          corresponding HMSL id
%   DrugLayout      cell array of matrix with doses, order corresponds to
%                       DrugNames; matrix size to plate size
%                       Can also be a 3D array with the third dimention
%                       being the DrugName dimension
%   nReps           Number of repeats (randomized replicates)
%
%
%   Optional arguments (name-value pairs):
%
%   Seed            initial seed, set to 0 for no randomization (default 1).
%                       For each repeat the seed is 'Seed+i-1'
%   stock_conc      in uM, default is 1e4 (= 10mM)
%   well_volume     in uL, default is 60 (for 384-well plates)
%   max_DMSOpct     maximum percent of DMSO dispensed (default 0.2%)
%   min_volume      step volume dispensed (in uL, default is 2e-5)
%   step_volume     step volume dispensed (in uL, default is 1.3e-5)
%   treated_wells   matrix of size plate_dims that display the
%                       treated wells; randomization will be within these
%                       wells (by default, use the convex rectangle of
%                       treatments defiend in DrugLayout)
%
%   Output: array of Design structure matching the DrugResponseSuite
%
%


%% inputs and parser
p = inputParser;

addParameter(p,'Seed',1,  @(x) isscalar(x) && isnumeric(x));
addParameter(p,'stock_conc',1e4, @(x) isvector(x) && isnumeric(x));     % in uM
addParameter(p,'well_volume',60, @(x) isscalar(x) && isnumeric(x));       % in uL
addParameter(p,'treated_wells',[], @(x) islogical(x));

% based on the specifications of the D300
addParameter(p,'min_volume', 1.3e-5, @(x) isscalar(x) && isnumeric(x)); % minimum volume of 13pl (in uL)
addParameter(p,'step_volume', 2e-5, @(x) isscalar(x) && isnumeric(x)); % minimum step of 20pl (in uL)
addParameter(p,'max_DMSOpct', .2, @(x) isscalar(x) && isnumeric(x)); % maximum 0.2% percent of DMSO

parse(p, varargin{:})
p = p.Results;
for i = {'Seed' 'stock_conc' 'treated_wells' ...
        'well_volume' 'min_volume' 'step_volume' 'max_DMSOpct'}
    eval([i{:} ' = p.' i{:} ';'])
end
% avoid too high level of DMSO (max is 0.2%)
max_volume = well_volume*max_DMSOpct/100;

%% checks
DrugNames = ToColumn(DrugNames);
assert((length(DrugNames)==length(HMSLids)) || isempty(HMSLids))
if iscell(DrugLayout)
    assert(length(DrugNames)==length(DrugLayout))
    plate_dims = unique(cell2mat(cellfun2(@(x) size(x)',DL))','rows');
    assert(size(plate_dims,1)==1)
else
    plate_dims = size(DrugLayout);
    if numel(plate_dims)==2, plate_dims = [plate_dims 1]; end
    assert(length(DrugNames)==plate_dims(3));
    plate_dims = plate_dims(1:2);
    DrugLayout = squeeze(num2cell(DrugLayout,[1 2]));
end
if length(plate_dims)==1
    plate_dims = sqrt(plate_dims*[1/1.5 1.5]);
end
assert(round(log2(plate_dims(1)))==log2(plate_dims(1)) && (plate_dims(1)*1.5==plate_dims(2)),...
    'Non standard plate size: %i x %i', plate_dims(1), plate_dims(2))


assert((length(DrugNames)==length(stock_conc)) || (length(stock_conc)==1))
if ~(all(stock_conc>=10 & stock_conc<=5e4))
    warning('Stock concentration should be in uM!')
    pause
end

if length(stock_conc)==1
    stock_conc = stock_conc*ones(length(DrugNames),1);
end
stock_conc = num2cell(ToColumn(stock_conc));


if Seed==0 && nReps>1
    warnprintf('Seed is set to 0 (i.e. no randomization) and Nreps>1')
    warnprintf('This means all replicates will have the same layout!')
    warnprintf('Confirm to proceed')
    pause
end


%% initialization
if ~isempty(treated_wells)
    assert(all(size(treated_wells)==plate_dims));
else
    allDrugs = reshape([DrugLayout{:}], [plate_dims length(DrugNames)]);
    treated_wells = any(allDrugs,3);
    [r,c] = find(treated_wells>0);
    treated_wells = false(plate_dims);
    treated_wells(min(r):max(r), min(c):max(c)) = true;
end

Drugs = struct('DrugName', DrugNames, 'HMSLid', ToColumn(HMSLids),...
    'stock_conc', stock_conc, 'layout', zeros(plate_dims));
Designs = struct('plate_dims', repmat({plate_dims}, nReps, 1), ...
    'treated_wells', repmat({treated_wells}, nReps, 1), ...
    'well_volume', repmat({well_volume}, nReps, 1), ...
    'Drugs', repmat({Drugs}, nReps, 1), 'Seed', num2cell(Seed-1+(1:nReps)'));

DrugLayout = cellfun2(@(x,y,z) round_Doses(x,y,z,'layout',min_volume, ...
    step_volume, max_volume, well_volume), DrugLayout, stock_conc, DrugNames);

Wellidx = find(treated_wells);
nWells = length(Wellidx);

%% assign the design for each replicate
for iR = 1:nReps

    s = RandStream('mt19937ar','Seed',Seed+iR-1);
    RandStream.setGlobalStream(s);
    if Seed>0
        idx = randperm(nWells); % randomize the order
    else
        idx = 1:nWells;
    end
    order = sortidx(idx','ascend'); % find the order for the treatment on the plate

    for iD = 1:length(DrugNames)
        Designs(iR).Drugs(iD).layout(Wellidx(order)) = ...
            DrugLayout{iD}(Wellidx);
    end


end


end
