function UpdatedClassifiedCells = CellCycleReClass(ClassifiedCells, ManualIdxs, SelectedWells, varargin)
% ClassifiedCells = CellCycleClassification(t_SingleCelldata, SingleCelldata, ...)
%
% inputs are :
%  t_SingleCelldata  -> metadata for the plates/wells
%  SingleCelldata    -> object-level data (array of structures with one field per
%                               channel, same order as t_SingleCelldata)
%
% parameter inputs are:
%  Channelnames -> name of the variables in SingleCelldata (structure with
%                       field LDR, DNA, and optional: EdU, pH3
%  savefolder   -> name to save images of the results
%%%   option to save only if a test fails ?? <<<<<<<<<<<<<<<----------------------
%  ConditionKeys-> name of keys to split conditions
%  PosCtrlLabel -> labels for the positive controls
%  NegCtrlLabel -> labels for the negative controls
%
% output is a structure with:
%  t_results       -> table with all calculated data and metadata
%  allCellIdentity -> cell cycle phase assigned to each cell
%                       (dead are -1; unclassified 0; G1 1; S 2; G2 3; M 4)
%  t_qc            -> table with summary of qc results by well
%  t_summary       -> table with summary of qc results by condition
%  plotResults     -> structure with the values and gates for plotting EdU/DNA
%
%% assign inputs and prepare tests

if islogical(ManualIdxs)
    assert(height(ClassifiedCells.t_results)==length(ManualIdxs))
    ManualIdxs = find(ManualIdxs);
    
elseif isnumeric(ManualIdxs) 
    assert(all(ManualIdxs>0 & ManualIdxs<=height(ClassifiedCells.t_results)))
    
else
    SelectedBarcodes = ManualIdxs;
    if ischar(SelectedBarcodes), SelectedBarcodes = {SelectedBarcodes}; end
    if ischar(SelectedWells), SelectedWells = {SelectedWells}; end
    if iscellstr(SelectedWells)
        SelectedWells = {SelectedWells};
    end
    assert(length(SelectedWells)==length(SelectedWells))
    SelectedWells
    ManualIdxs = false(height(ClassifiedCells.t_results),1);
    for iB = 1:length(SelectedBarcodes)
        ManualIdxs(ClassifiedCells.t_results.Barcode==SelectedBarcodes{iB} & ...
            ismember(ClassifiedCells.t_results.Well,SelectedWells{iB})) = true;
        
    end
    ManualIdxs = find(ManualIdxs);
end

p = parseCellCycleInputs(varargin{:});

% flag for the modular approach
useEdU = ~isempty(ClassifiedCells.plotResults(1).EdU);
usepH3 = ~isempty(ClassifiedCells.plotResults(1).pH3);

% folder for storing data
if ~isempty(p.savefolder)
    grp_savefolder = [p.savefolder '/manual/'];
    mkdir(grp_savefolder)
else
    savefig = ''; % this will bypass saving
end

% unwrap the structure
t_results = ClassifiedCells.t_results;
t_qc = ClassifiedCells.t_qc;
allCellIdentity = ClassifiedCells.allCellIdentity;
plotResults = ClassifiedCells.plotResults;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through the selected conditions

fprintf('\nReprocessing %i wells: ', length(ManualIdxs));
for iW = 1:length(ManualIdxs)
    
    if ismember(t_results.pert_type(ManualIdxs(iW)), p.NegCtrlLabel)
        welltag = 'Negctl_';
    elseif ismember(t_results.pert_type(ManualIdxs(iW)), p.PosCtrlLabel)
        welltag = 'Posctl_';
    else
        welltag = 'Trt_';
    end
    fprintf('\n  % 2i/%i: %s - %s: %s, c=%.2g (%s)', iW, length(ManualIdxs), ...
        t_results.Barcode(ManualIdxs(iW)), t_results.Well(ManualIdxs(iW)), ...
        t_results.DrugName(ManualIdxs(iW)), t_results.Conc(ManualIdxs(iW)), welltag(1:end-1));
    logtxt = 'MANUAL:';
    
    
    LDR = ClassifiedCells.plotResults(ManualIdxs(iW)).LDR;
    DNA = ClassifiedCells.plotResults(ManualIdxs(iW)).DNA;
    
    if isempty(LDR), continue, end % display warning or print in log? <-----------------
    
    if ~isempty(p.savefolder), savefig = [grp_savefolder ...
            welltag char(t_results.Well(ManualIdxs(iW))) '_LDR_manual.jpg'];
    end
    [LiveCells, DeadCells, plotResults(ManualIdxs(iW)).LDRGates, ...
        plotResults(ManualIdxs(iW)).DNAGates, CellOutcome, ~, ~, ltxt] = DeadCellFilter(LDR, DNA, ...
        'savefig', savefig, 'interactive', true);
    % store the results
    logtxt = [logtxt ' ' ltxt];
    allCellIdentity{ManualIdxs(iW)} = CellOutcome;
    t_results(ManualIdxs(iW), {'LiveCells' 'DeadCells'}) = {LiveCells, DeadCells};
    % check and report outcome (special for toxic positive control)
    t_qc.PassFracDead(ManualIdxs(iW)) = true;
    
    if useEdU
        EdU = ClassifiedCells.plotResults(ManualIdxs(iW)).EdU;
        if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                welltag char(t_results.Well(ManualIdxs(iW))) '_EdU_manual.jpg'];
        end
        [CCPks, CCfrac, ...
            plotResults(ManualIdxs(iW)).DNAGates, plotResults(ManualIdxs(iW)).EdUGates, ...
            CellIdentity, ...
            plotResults(ManualIdxs(iW)).logDNA, plotResults(ManualIdxs(iW)).logEdU,~,~,ltxt] = ...
            CCphases(DNA(CellOutcome==1), EdU(CellOutcome==1), ...
            'savefig', savefig, 'interactive', true);
        logtxt = [logtxt '; ' ltxt];
        
        if usepH3
            pH3 = ClassifiedCells.plotResults(ManualIdxs(iW)).pH3;
            if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                    welltag char(t_results.Well(ManualIdxs(iW))) '_pH3_manual.jpg'];
            end
            [CCfrac, CellIdentity, plotResults(ManualIdxs(iW)).H3cutoff, ~, ~, ltxt] = ...
                pH3Filter(pH3(CellOutcome==1), ...
                CellIdentity, 'savefig', savefig, 'interactive', true);
            logtxt = [logtxt '; ' ltxt];
        end
        
        % store the results
        allCellIdentity{ManualIdxs(iW)}(CellOutcome==1) = CellIdentity;
        t_results(ManualIdxs(iW), {'CCfrac' 'CCPks'}) = {CCfrac {CCPks}};
        % check consistency and report outcome
        PassUnclass = mean(CellIdentity==0) < p.TestCutoffs.Unclass;
        PassCCphase = true;
        PassPksConsist = true;
        
        t_qc(ManualIdxs(iW), ...
            {'PassUnclass' 'PassCCphase' 'PassPksConsist'}) = ...
            {PassUnclass PassCCphase PassPksConsist};
        t_qc.log{ManualIdxs(iW)} = logtxt;
    end
    
end
fprintf(' done\n')



%%
t_summary = collapse(t_qc, @mean, ...
    'keyvars', unique([{'Barcode', 'CellLine'} p.ConditionKeys 'pert_type'],'stable'), ...
    'valvars', intersect(varnames(t_qc), {'PassFracDead' 'PassDeadConsist' 'PassUnclass' ...
    'PassCCphase' 'PassPksConsist' 'PassCCConsist'}, 'stable'));


UpdatedClassifiedCells = struct('t_results', {t_results}, 'allCellIdentity', {allCellIdentity}, ...
    't_qc', {t_qc}, 't_summary', {t_summary}, 'plotResults', {plotResults});