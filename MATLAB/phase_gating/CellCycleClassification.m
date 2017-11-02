function ClassifiedCells = CellCycleClassification(t_SingleCelldata, SingleCelldata, varargin)
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
% output is a structure ClassifiedCells with fields:
%  t_results       -> table with all calculated data and metadata
%  allCellIdentity -> cell cycle phase assigned to each cell 
%                       (dead are -1; unclassified 0; G1 1; S 2; G2 3; M 4)
%  t_qc            -> table with summary of qc results by well
%  t_summary       -> table with summary of qc results by condition
%  plotResults     -> structure with the values and gates for plotting EdU/DNA
%
%% assign inputs and prepare tests

assert(height(t_SingleCelldata)==length(SingleCelldata))

p = parseCellCycleInputs(varargin{:});


% flag for the modular approach
useEdU = ~isempty(p.Channelnames.EdU);
usepH3 = ~isempty(p.Channelnames.pH3);

% save the results
plotResults = struct('LDR', {SingleCelldata.(p.Channelnames.LDR)}', ...
    'LDRGates', repmat({[]}, height(t_SingleCelldata),1), ...
    'DNApreGates', repmat({[]}, height(t_SingleCelldata),1), ...
    'DNA', {SingleCelldata.(p.Channelnames.DNA)}', ...
    'EdU', {SingleCelldata.(p.Channelnames.EdU)}', ...
    'logDNA', repmat({[]}, height(t_SingleCelldata),1), ...
    'logEdU', repmat({[]}, height(t_SingleCelldata),1), ...
    'DNAGates', repmat({[]}, height(t_SingleCelldata),1), ...
    'EdUGates', repmat({[]}, height(t_SingleCelldata),1), ...
    'pH3', {SingleCelldata.(p.Channelnames.pH3)}', ...
    'pH3cutoff', repmat({[]}, height(t_SingleCelldata),1));
    
%% %%%%%%%%%%%%%%%%%%
% start the analysis

% get all different conditions
default_metadata = {'CellLine' 'DrugName' 'Conc' 'DrugName2' 'Conc2' 'Time'};
t_groups = unique(t_SingleCelldata(:,unique([{'Barcode'} p.ConditionKeys],'stable')));

% start the result report
t_results = [t_SingleCelldata(:,unique([ ...
    intersect(default_metadata, varnames(t_SingleCelldata), 'stable'), ...
    p.ConditionKeys {'pert_type' 'Barcode' 'Well'}],'stable')) ...
    array2table(NaN(height(t_SingleCelldata),2), 'variablenames', ...
    {'LiveCells' 'DeadCells'})];
t_qc = [t_results(:, unique([{'Barcode', 'CellLine'} p.ConditionKeys 'pert_type'],'stable')) ...
    array2table(false(height(t_SingleCelldata),2),'variablenames', ...
    {'PassFracDead' 'PassDeadConsist'}) ...
    table(repmat({''},height(t_SingleCelldata),1),'variablenames', {'log'})];

if useEdU
    t_results = [t_results ...
        table(NaN(height(t_SingleCelldata),4+usepH3), ...
        repmat({NaN(3,2)}, height(t_SingleCelldata),1), ...
        'variablenames', {'CCfrac' 'CCPks'})];
    t_qc = [t_qc ...
        array2table(false(height(t_SingleCelldata),4),'variablenames', ...
        {'PassUnclass' 'PassCCphase' 'PassCCConsist' 'PassPksConsist'})];
end

% structure to store all results per well
allCellIdentity = cell(height(t_results),1);
CellIdentity = [];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through the different conditions

for iGr = 1:height(t_groups)
    fprintf('\n%i of %i: %s\n-----------------------------------\n', ...
        iGr, height(t_groups), strjoin(table2cellstr(t_groups(iGr,:),0),' '));
    
    % folder for storing data
    if ~isempty(p.savefolder)
        grp_savefolder = [p.savefolder '/' ...
            strjoin(table2cellstr(t_groups(iGr,:),0),'_') '/'];
        mkdir(grp_savefolder)
    else
        savefig = ''; % this will bypass plotting and saving
    end
    
    %% %%%%%% all controls first %%%%%%%%%%%
    fprintf('\tAll controls\n');
    % first get all controls to define the dead cell gating
    NegCtrlidx = find(eqtable(t_groups(iGr,:), t_SingleCelldata) & ...
        ismember(t_SingleCelldata.pert_type, p.NegCtrlLabel));
    
    PosCtrlidx = find(eqtable(t_groups(iGr,:), t_SingleCelldata) & ...
        ismember(t_SingleCelldata.pert_type, p.PosCtrlLabel));
    
    % if no positive control, use all wells            <<<<<<<<< -- not sure if right --------
    if any(PosCtrlidx), AllCtrlidx = [NegCtrlidx; PosCtrlidx]; else
        AllCtrlidx = find(eqtable(t_groups(iGr,:), t_SingleCelldata)); end
    %  <<<<<<<<< -- maybe safer to only use the negative controls --------
    
    LDR = vertcat(SingleCelldata(AllCtrlidx).(p.Channelnames.LDR));
    DNA = vertcat(SingleCelldata(AllCtrlidx).(p.Channelnames.DNA));
    
    if ~isempty(p.savefolder), savefig = [grp_savefolder 'All_ctl_LDR.jpg']; end
    [~, ~, ~, ~, ~, LDRlims, DNAlims] = DeadCellFilter(LDR, DNA, ...
        'savefig', savefig);
    
    
    % run on all negative controls
    LDR = vertcat(SingleCelldata(NegCtrlidx).(p.Channelnames.LDR));
    DNA = vertcat(SingleCelldata(NegCtrlidx).(p.Channelnames.DNA));
    if ~isempty(p.savefolder), savefig = [grp_savefolder 'All_Negctl_LDR.jpg']; end
    [~, RefFracDead, ~, ~, CellOutcome] = DeadCellFilter(LDR, DNA, ...
        'savefig', savefig, 'LDRlims', LDRlims, 'DNAlims', DNAlims);
    RefFracDead = RefFracDead/length(LDR);
    meanNcells = sum(CellOutcome==1)/length(NegCtrlidx);
    if useEdU
        EdU = vertcat(SingleCelldata(NegCtrlidx).(p.Channelnames.EdU));
        if ~isempty(p.savefolder), savefig = [grp_savefolder 'All_Negctl_EdU.jpg']; end
        [RefCCPks, RefCCfrac,~,~,CellIdentity,~,~,~, EdUlims] = CCphases(DNA(CellOutcome==1), ...
            EdU(CellOutcome==1), 'savefig', savefig, 'DNAlims', DNAlims);
    end
    if usepH3
        pH3 = vertcat(SingleCelldata(NegCtrlidx).(p.Channelnames.pH3));
        if ~isempty(p.savefolder), savefig = [grp_savefolder 'All_Negctl_pH3.jpg']; end
        [RefCCfrac, ~, RefpH3cutoff, pH3lims] = pH3Filter(pH3(CellOutcome==1), CellIdentity, ...
            'savefig', savefig);
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check how it applies for each negative control
    fprintf('\tNeg controls (%i): ', length(NegCtrlidx));
    for iW = 1:length(NegCtrlidx)
        fprintf(' %s', char(t_SingleCelldata.Well(NegCtrlidx(iW))));
        
        LDR = SingleCelldata(NegCtrlidx(iW)).(p.Channelnames.LDR);
        DNA = SingleCelldata(NegCtrlidx(iW)).(p.Channelnames.DNA);
        
        if isempty(LDR), continue, end % display warning or print in log? <-----------------
        
        if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                'Negctl_' char(t_SingleCelldata.Well(NegCtrlidx(iW))) '_LDR.jpg'];
        end
        [LiveCells, DeadCells, plotResults(NegCtrlidx(iW)).LDRGates, ...
            plotResults(NegCtrlidx(iW)).DNAGates, CellOutcome,~,~,ltxt] = ...
            DeadCellFilter(LDR, DNA, ...
            'LDRlims', LDRlims, 'DNAlims', DNAlims, ...
            'savefig', savefig);
        logtxt = ['Auto: ' ltxt];
        % store the results
        allCellIdentity{NegCtrlidx(iW)} = CellOutcome;
        t_results(NegCtrlidx(iW), {'LiveCells' 'DeadCells'}) = {LiveCells, DeadCells};
        % check consistency and report outcome
        PassFracDead = DeadCells/length(LDR) < p.TestCutoffs.FracDead;
        PassDeadConsist = abs(DeadCells/length(LDR) - RefFracDead) < p.TestCutoffs.DeadConsist;
        
        t_qc(NegCtrlidx(iW), {'PassFracDead' 'PassDeadConsist'}) = ...
            {PassFracDead PassDeadConsist};
        
        if useEdU
            EdU = SingleCelldata(NegCtrlidx(iW)).(p.Channelnames.EdU);
            if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                    'Negctl_' char(t_SingleCelldata.Well(NegCtrlidx(iW))) '_EdU.jpg'];
            end
            [CCPks, CCfrac, ...
                plotResults(NegCtrlidx(iW)).DNAGates, plotResults(NegCtrlidx(iW)).EdUGates, ...
                CellIdentity, ...
                plotResults(NegCtrlidx(iW)).logDNA, plotResults(NegCtrlidx(iW)).logEdU,~,~,ltxt] = ...
                CCphases(DNA(CellOutcome==1), EdU(CellOutcome==1), ...
                'savefig', savefig, 'DNAlims', DNAlims, 'EdUlims', EdUlims);
            logtxt = [logtxt '; ' ltxt];
            
            if usepH3
                pH3 = SingleCelldata(NegCtrlidx(iW)).(p.Channelnames.pH3);
                if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                        'Negctl_' char(t_SingleCelldata.Well(NegCtrlidx(iW))) '_pH3.jpg'];
                end
                [CCfrac, CellIdentity, plotResults(NegCtrlidx(iW)).H3cutoff, ~, ...
                    plotResults(NegCtrlidx(iW)).logpH3, ltxt] = pH3Filter(pH3(CellOutcome==1), ...
                    CellIdentity, 'savefig', savefig, 'pH3lims', pH3lims, ...
                    'pH3cutoff', RefpH3cutoff);
                logtxt = [logtxt '; ' ltxt];
            end
            
            % store the results
            allCellIdentity{NegCtrlidx(iW)}(CellOutcome==1) = CellIdentity;
            t_results(NegCtrlidx(iW), {'CCfrac' 'CCPks'}) = {CCfrac {CCPks}};
            % check consistency and report outcome
            PassUnclass = mean(CellIdentity==0) < p.TestCutoffs.Unclass;
            PassCCphase = all(CCfrac(1:3)>.05) & CCfrac(end-1) > p.TestCutoffs.CCphase;
            PassCCConsist = max(abs(CCfrac - RefCCfrac)) < p.TestCutoffs.CCConsist;
            PassPksConsist = max(abs(CCPks(:,1)-RefCCPks(:,1))) < p.TestCutoffs.PksConsist;
            
            t_qc(NegCtrlidx(iW), ...
                {'PassUnclass' 'PassCCphase' 'PassCCConsist' 'PassPksConsist'}) = ...
                {PassUnclass PassCCphase PassCCConsist PassPksConsist};
            t_qc.log{NegCtrlidx(iW)} = logtxt;
        end
    end
    fprintf(' done\n')
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check how it applies for each positive control
    
    fprintf('\tPos controls (%i): ', length(PosCtrlidx));
    for iW = 1:length(PosCtrlidx)
        fprintf(' %s', char(t_SingleCelldata.Well(PosCtrlidx(iW))));
        
        LDR = SingleCelldata(PosCtrlidx(iW)).(p.Channelnames.LDR);
        DNA = SingleCelldata(PosCtrlidx(iW)).(p.Channelnames.DNA);
        plotResults(PosCtrlidx(iW)).LDR = LDR;
        plotResults(PosCtrlidx(iW)).DNA = DNA;
        
        if isempty(LDR), continue, end % display warning or print in log? <-----------------
        
        if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                'Posctl_' char(t_SingleCelldata.Well(PosCtrlidx(iW))) '_LDR.jpg'];
        end
        [LiveCells, DeadCells, plotResults(PosCtrlidx(iW)).LDRGates, ...
            plotResults(PosCtrlidx(iW)).DNAGates, CellOutcome,~,~,ltxt] = ...
            DeadCellFilter(LDR, DNA, ...
            'LDRlims', LDRlims, 'DNAlims', DNAlims, ...
            'savefig', savefig);
        logtxt = ['Auto: ' ltxt];
        % store the results
        allCellIdentity{PosCtrlidx(iW)} = CellOutcome;
        t_results(PosCtrlidx(iW), {'LiveCells' 'DeadCells'}) = {LiveCells, DeadCells};
        % check and report outcome (special for toxic positive control)
        t_qc.PassFracDead(PosCtrlidx(iW)) = (DeadCells/length(LDR)) > ...
            (RefFracDead + .5*(t_SingleCelldata.pert_type(PosCtrlidx(iW))=='ctl_toxic'));
        
        if useEdU
            EdU = SingleCelldata(PosCtrlidx(iW)).(p.Channelnames.EdU);
            plotResults(PosCtrlidx(iW)).EdU = EdU;
            if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                    'Posctl_' char(t_SingleCelldata.Well(PosCtrlidx(iW))) '_EdU.jpg'];
            end
            [CCPks, CCfrac, ...
                plotResults(PosCtrlidx(iW)).DNAGates, plotResults(PosCtrlidx(iW)).EdUGates, ...
                CellIdentity, ...
                plotResults(PosCtrlidx(iW)).logDNA, plotResults(PosCtrlidx(iW)).logEdU,~,~,ltxt] = ...
                CCphases(DNA(CellOutcome==1), EdU(CellOutcome==1), ...
                'savefig', savefig, 'DNAlims', DNAlims, 'EdUlims', EdUlims);
            logtxt = [logtxt '; ' ltxt];
            
            if usepH3
                pH3 = SingleCelldata(PosCtrlidx(iW)).(p.Channelnames.pH3);
                if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                        'Posctl_' char(t_SingleCelldata.Well(PosCtrlidx(iW))) '_pH3.jpg'];
                end
                if sum(CellOutcome==1)>meanNcells*p.CellFrac_Adpative_pH3cutoff %% check if low cell count
                    temp_pH3cutoff = RefpH3cutoff;
                else
                    temp_pH3cutoff = [];
                end
                [CCfrac, CellIdentity, plotResults(PosCtrlidx(iW)).H3cutoff, ~, ...
                    plotResults(PosCtrlidx(iW)).logpH3, ltxt] = pH3Filter(pH3(CellOutcome==1), ...
                    CellIdentity, 'savefig', savefig, 'pH3lims', pH3lims, ...
                    'pH3cutoff', temp_pH3cutoff);
                logtxt = [logtxt '; ' ltxt];
                    
            end
            
            % store the results
            allCellIdentity{PosCtrlidx(iW)}(CellOutcome==1) = CellIdentity;
            t_results(PosCtrlidx(iW), {'CCfrac' 'CCPks'}) = {CCfrac {CCPks}};
            % check consistency and report outcome
            PassUnclass = mean(CellIdentity==0) < .05;
            switch t_SingleCelldata.pert_type(PosCtrlidx(iW))
                case 'ctl_G1'
                    PassCCphase = CCfrac(1) > p.TestCutoffs.CCphase_pos(RefCCfrac(1));
                case 'ctl_S'
                    PassCCphase = CCfrac(2) > p.TestCutoffs.CCphase_pos(RefCCfrac(2));
                case 'ctl_G2'
                    PassCCphase = CCfrac(3) > p.TestCutoffs.CCphase_pos(RefCCfrac(3));                    
                case 'ctl_M'
                    if usepH3
                        PassCCphase = CCfrac(4) > p.TestCutoffs.CCphase_pos(RefCCfrac(4));
                    else, PassCCphase = true;  end
                otherwise
                    PassCCphase = true;  % <<<<<<<<<<<<------------ need to check
            end
            PassPksConsist = max(abs(CCPks(:,1)-RefCCPks(:,1))) < p.TestCutoffs.PksConsist;
            
            t_qc(PosCtrlidx(iW), ...
                {'PassUnclass' 'PassCCphase' 'PassPksConsist'}) = ...
                {PassUnclass PassCCphase PassPksConsist};
            t_qc.log{PosCtrlidx(iW)} = logtxt;
        end
    end
    fprintf(' done\n')
    
    %% check consistency among the positive controls with the same treatments
    t_Posctrl = unique(t_SingleCelldata(PosCtrlidx, ...
        unique([{'Barcode', 'CellLine'} p.ConditionKeys {'DrugName' 'Conc'}])));
    for iP=1:height(t_Posctrl)
        idx = find(eqtable(t_Posctrl(iP,:), t_SingleCelldata));
        FracDead = t_results.DeadCells(idx)./ ...
            (t_results.DeadCells(idx) + t_results.LiveCells(idx));
        t_qc.PassDeadConsist(idx) = abs(FracDead - mean(FracDead)) < p.TestCutoffs.DeadConsist;
        
        CCfrac = t_results.CCfrac(idx,:);
        t_qc.PassCCConsist(idx) = max(abs(CCfrac - ...
            repmat(mean(CCfrac,1),length(idx),1)),[],2) < 2*p.TestCutoffs.CCConsist;
    end
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % work out all the other wells in the same condition
    Trtidx = find(eqtable(t_groups(iGr,:), t_SingleCelldata) & ...
        ~ismember(t_SingleCelldata.pert_type, [p.NegCtrlLabel p.PosCtrlLabel]));
    fprintf('\tTreatment wells (%i): ', length(Trtidx));
    for iW = 1:length(Trtidx)
        if mod(iW, ceil(length(Trtidx)/8))==0, fprintf(' %i', iW); end
        
        LDR = SingleCelldata(Trtidx(iW)).(p.Channelnames.LDR);
        DNA = SingleCelldata(Trtidx(iW)).(p.Channelnames.DNA);
        plotResults(Trtidx(iW)).LDR = LDR;
        plotResults(Trtidx(iW)).DNA = DNA;
        
        if isempty(LDR), continue, end % display warning or print in log? <-----------------
        
        if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                'Trt_' char(t_SingleCelldata.Well(Trtidx(iW))) '_LDR.jpg'];
        end
        [LiveCells, DeadCells, plotResults(Trtidx(iW)).LDRGates, plotResults(Trtidx(iW)).DNAGates, CellOutcome,~,~,ltxt] = ...
            DeadCellFilter(LDR, DNA, 'LDRlims', LDRlims, 'DNAlims', DNAlims, ...
            'savefig', savefig);
        logtxt = ['Auto: ' ltxt];
        
        % store the results
        allCellIdentity{Trtidx(iW)} = CellOutcome;
        t_results(Trtidx(iW), {'LiveCells' 'DeadCells'}) = {LiveCells, DeadCells};
        % check and report outcome (special for toxic positive control)
        t_qc.PassFracDead(Trtidx(iW)) = true;
        
        if useEdU
            EdU = SingleCelldata(Trtidx(iW)).(p.Channelnames.EdU);
            plotResults(Trtidx(iW)).EdU = EdU;
            if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                    'Trt_' char(t_SingleCelldata.Well(Trtidx(iW))) '_EdU.jpg'];
            end
            [CCPks, CCfrac, ...
                plotResults(Trtidx(iW)).DNAGates, plotResults(Trtidx(iW)).EdUGates, ...
                CellIdentity, ...
                plotResults(Trtidx(iW)).logDNA, plotResults(Trtidx(iW)).logEdU,~,~,ltxt] = ...
                CCphases(DNA(CellOutcome==1), EdU(CellOutcome==1), ...
                'savefig', savefig, 'DNAlims', DNAlims, 'EdUlims', EdUlims);
            logtxt = [logtxt '; ' ltxt];
            
            if usepH3
                pH3 = SingleCelldata(Trtidx(iW)).(p.Channelnames.pH3);
                logResults(Trtidx(iW)).pH3 = pH3;
                if ~isempty(p.savefolder), savefig = [grp_savefolder ...
                        'Trt_' char(t_SingleCelldata.Well(Trtidx(iW))) '_pH3.jpg'];
                end
                if sum(CellOutcome==1)>meanNcells*p.CellFrac_Adpative_pH3cutoff %% check if low cell count
                    temp_pH3cutoff = RefpH3cutoff;
                else
                    temp_pH3cutoff = [];
                end
                [CCfrac, CellIdentity, plotResults(Trtidx(iW)).H3cutoff, ~, ...
                    ~, ltxt] = pH3Filter(pH3(CellOutcome==1), ...
                    CellIdentity, 'savefig', savefig, 'pH3lims', pH3lims, ...
                    'pH3cutoff', temp_pH3cutoff);
                logtxt = [logtxt '; ' ltxt];
            end
            
            % store the results
            allCellIdentity{Trtidx(iW)}(CellOutcome==1) = CellIdentity;
            t_results(Trtidx(iW), {'CCfrac' 'CCPks'}) = {CCfrac {CCPks}};
            % check consistency and report outcome
            PassUnclass = mean(CellIdentity==0) < p.TestCutoffs.Unclass;
            PassCCphase = true;  
            PassPksConsist = max(abs(CCPks(:,1)-RefCCPks(:,1))) < p.TestCutoffs.PksConsist;
            
            t_qc(Trtidx(iW), ...
                {'PassUnclass' 'PassCCphase' 'PassPksConsist'}) = ...
                {PassUnclass PassCCphase PassPksConsist};
            t_qc.log{Trtidx(iW)} = logtxt;
        end
    end
    
    %% check consistency among the wells with same treatments
    t_trt = unique(t_SingleCelldata(Trtidx,  ...
        unique([{'Barcode', 'CellLine'} p.ConditionKeys {'DrugName' 'Conc'}])));
    for iP = 1:height(t_trt)
        idx = find(eqtable(t_trt(iP,:), t_SingleCelldata));
        FracDead = t_results.DeadCells(idx)./ ...
            (t_results.DeadCells(idx) + t_results.LiveCells(idx));
        t_qc.PassDeadConsist(idx) = abs(FracDead - nanmean(FracDead)) < p.TestCutoffs.DeadConsist;
        
        CCfrac = t_results.CCfrac(idx,:);
        t_qc.PassCCConsist(idx) = max(abs(CCfrac - ...
            repmat(nanmean(CCfrac,1),length(idx),1)),[],2) < 2*p.TestCutoffs.CCConsist;
    end
    
    fprintf('\tdone\n')
end


%%
t_summary = collapse(t_qc, @mean, ...
    'keyvars', unique([{'Barcode', 'CellLine'} p.ConditionKeys 'pert_type'],'stable'), ...
    'valvars', intersect(varnames(t_qc), {'PassFracDead' 'PassDeadConsist' 'PassUnclass' ...
    'PassCCphase' 'PassPksConsist' 'PassCCConsist'}, 'stable'));


ClassifiedCells = struct('t_results', {t_results}, 'allCellIdentity', {allCellIdentity}, ...
    't_qc', {t_qc}, 't_summary', {t_summary}, 'plotResults', {plotResults});