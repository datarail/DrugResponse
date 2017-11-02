function t_annotated = Annotate_CellCountData(t_data, folder, fields)
% t_annotated = Annotate_CellCountData(t_data, folder, fields)
%
%   annotate the data using the treatment files
%   adding the fields:
%       - DrugName
%       - Conc
%       - selected additional fields as given in input varaible 'fields' (by
%           default ALL 'Perturbations' in the array of structures Design)
%
%   variable folder is to specify where the TreatmentFiles are stored.
%
fprintf('Annotate Cell count data with plate designs:\n');

if ~exist('folder','var') || isempty(folder)
    folder = '';
end

%% import all treatment files and read the designs
Trtfiles = setdiff(cellstr(unique(t_data.TreatmentFile)),'-');

Designs = cell(length(Trtfiles),1);
DesignNumbers = cell(length(Trtfiles),1);
correct_barcode = Designs;
allDesigns = {};
for iTf = 1:length(Trtfiles)
    assert(exist(fullfile(folder, Trtfiles{iTf}), 'file')>0, ...
        'Treatment file %s missing in folder %s', Trtfiles{iTf}, folder)
    [~,~,ext] = fileparts(Trtfiles{iTf});
    fprintf('\tLoading %s\n', Trtfiles{iTf})

    DesignNumbers{iTf} = unique(t_data.DesignNumber(t_data.TreatmentFile==Trtfiles{iTf}));
    if strcmp(ext,'.mat')
        temp = load(fullfile(folder, Trtfiles{iTf}));
        try Designs{iTf} = temp.Design; catch, Designs{iTf} = temp.Designs; end
    elseif strcmp(ext,'.hpdd')
        [Designs{iTf}, correct_barcode{iTf}] = hpdd_importer(fullfile(folder, Trtfiles{iTf}));
        % because of redundant plates, barcodes have to be reassigned
        % barcode found in the hpdd files will overwrite other barcode
    else
        try % assume a tsv file
            Designs{iTf} = TextDesignFile_importer(fullfile(folder, Trtfiles{iTf}));
        catch err
            warnprintf(['expecting a .mat, .hpdd or a formated ' ...
                'tab-separated file as TreatmentFile\n'])
            rethrow(err)
        end
    end
    allDesigns = [allDesigns;
        num2cell(Designs{iTf}(DesignNumbers{iTf}))];
end

%% look up all the drugs and perturbations in the designs

Ndrugs = 1;
Perturbations = {};
DrugNames = {};
t_HMSLids = table;
for iD=1:numel(allDesigns)

    % check for multiple drugs in the same well
    if ~isfield(allDesigns{iD},'Drugs') || isempty(allDesigns{iD}.Drugs)
        continue
    end

    DrugConc = reshape([allDesigns{iD}.Drugs.layout], [allDesigns{iD}.plate_dims ...
        length(allDesigns{iD}.Drugs)]);
    DrugNames = unique([DrugNames {allDesigns{iD}.Drugs.DrugName}],'stable');
    if isfield(allDesigns{iD}.Drugs,'HMSLid')
        t_HMSLids = [ t_HMSLids; table({allDesigns{iD}.Drugs.DrugName}', ...
            {allDesigns{iD}.Drugs.HMSLid}','VariableNames', {'DrugName' 'HMSLid'})];
    end
    if any(any(sum(DrugConc>0,3)>Ndrugs))
        Ndrugs = max(max(sum(DrugConc>0,3)));
        warnprintf('some wells have %i drugs, additional columns in output', ...
            Ndrugs)
    end

    % store all possible perturbations
    if isfield(allDesigns{iD}, 'Perturbations')
        Perturbations = unique([Perturbations {allDesigns{iD}.Perturbations.Name}], 'stable');
    end
end
t_HMSLids = unique(t_HMSLids);


%% declare the variables


% this is not the optimal way of storing multiple drugs because of the
% hierarcy between DrugName and Conc as well as the redudancy and possible
% swapping between Drug1 and Drug2 ; it makes matching between condition hard


if exist('fields','var')
    assert(all(ismember(fields, Perturbations)), ...
        'Not all ''fields'' found as perturbations in the design files')
else
    fields = Perturbations;
end
datafields = cell(1, length(fields));


%%


t_annotated = table;

for iTf = 1:length(Trtfiles)

    DesignNumbers = unique(t_data.DesignNumber(t_data.TreatmentFile==Trtfiles{iTf}));
    fprintf('\tDesign %s:\n', Trtfiles{iTf} )
    for iDN = 1:length(DesignNumbers)
        if ~isempty(correct_barcode{iTf})
            DNidx = correct_barcode{iTf}.DesignNumber(DesignNumbers(iDN));
            assert(all(t_data.Barcode(t_data.DesignNumber==DesignNumbers(iDN))== ...
                correct_barcode{iTf}.Barcode(DesignNumbers(iDN))), ...
                'Mismatch between barcodes in the hpdd file and the DesignNumber')
            fprintf('\thpdd exp %i -> design %i\n', DesignNumbers(iDN), DNidx);

        else
            DNidx = DesignNumbers(iDN);
        end
        t_design = DrugDesignToTable(Designs{iTf}(DNidx), fields, DrugNames);
        idx = find(t_data.TreatmentFile==Trtfiles{iTf} & t_data.DesignNumber==DesignNumbers(iDN));

        % could be extended to any variable that is found in both Plate
        % Info and Design
        Common_vars = setdiff(intersect(varnames(t_data), varnames(t_design)), 'Well');
        for iV = 1:length(Common_vars)
%         if isvariable(t_design, 'CellLine')
            warnprintf('%s is a perurbation; replacing what is found in the PlateInfo: %s', ...
                Common_vars{iV}, strjoin(cellstr(unique(t_data.CellLine(idx))'),';'))
        end
        varToKeep = setdiff(varnames(t_data), Common_vars);
%             varToKeep = setdiff(varnames(t_data), 'CellLine');
%         else
%             varToKeep = varnames(t_data);
%         end

        [temp, ia] = innerjoin(t_data(idx,varToKeep), t_design, 'keys', 'Well');
        if height(temp)<length(idx)
            warnprintf('Some wells (%s) have no annotations in file %s --> ignored', ...
                strjoin(cellstr(unique(t_data.Well(idx(setdiff(1:length(idx),ia))))'),', '), ...
                Trtfiles{iTf})
        end
        temp = AddDrugNameColumn(temp, Ndrugs);
        temp.Untrt(temp.pert_type=='Untrt') = true;

        if ~isempty(t_annotated)
            % new vars found in temp only
            newvars = setdiff(varnames(temp), varnames(t_annotated));
            for i=1:length(newvars)
                if isnumeric(temp.(newvars{i}))
                    t_annotated = [t_annotated ...
                        table(NaN(height(t_annotated),1), 'variablenames', newvars(i))];
                else
                    t_annotated = [t_annotated ...
                        table(repmat({'-'},height(t_annotated),1), 'variablenames', ...
                        newvars(i))];
                end
            end
            %  vars in annotated not found in temp
            newvars = setdiff(varnames(t_annotated), varnames(temp));
            for i=1:length(newvars)
                if isnumeric(t_annotated.(newvars{i}))
                    temp = [temp ...
                        table(NaN(height(temp),1), 'variablenames', newvars(i))];
                else
                    temp = [temp ...
                        table(repmat({'-'},height(temp),1), 'variablenames', ...
                        newvars(i))];
                end
            end
                
            t_annotated = [t_annotated; temp(:,varnames(t_annotated))];
        else
            t_annotated = temp;
        end

    end

end

% fill up the columns for the untreated plates before merging
Untrtidx = t_data.TreatmentFile=='-';
if any(Untrtidx)
    NDrugs = sum(strfindcell(varnames(t_annotated),'DrugName')==1);
    pert_type = repmat({'Untrt'}, sum(Untrtidx),1);

    if isempty(t_annotated)
        othervars = varnames(t_data);
        newvars = {};
    else
        newvars = setdiff(varnames(t_annotated), [varnames(t_data) {'pert_type'}], 'stable');
        othervars = intersect(varnames(t_annotated), varnames(t_data), 'stable');
    end
    temp = table;
    for iD=1:NDrugs
        % add the HMSLid by default; it will be filtered afterwards
        temp = [temp cell2table(repmat({'-'}, sum(Untrtidx),2), 'VariableName', ...
            {sprintf('DrugName%i', iD(iD>1)), sprintf('HMSLid%i', iD(iD>1))}) ...
            table(zeros(sum(Untrtidx),1), 'VariableName', {sprintf('Conc%i', iD(iD>1))})];
    end

    addvars = setdiff(newvars, varnames(temp));
    for i=1:length(addvars)
        %%% need to treat the case of SeedingNumber
        if strcmp(newvars{i}, 'SeedingNumber')
            error('Broadcasting of seeding number not implemented')
        end
        if isnumeric(t_annotated.(addvars{i}))
            temp = [temp table(zeros(sum(Untrtidx),1), 'VariableName', addvars(i))];
        else
            temp = [temp cell2table(repmat({'-'}, sum(Untrtidx),1), ...
                'VariableName', addvars(i))];
        end
    end

    if isempty(t_annotated)
        t_annotated = [t_data(Untrtidx,othervars) temp(:,newvars) table(pert_type)];
    else
        t_annotated = [
            [t_data(Untrtidx,othervars) temp(:,newvars) table(pert_type)]
            t_annotated(:,othervars) t_annotated(:,newvars)  t_annotated(:,'pert_type')];
    end
end

warnassert(height(t_annotated)==height(t_data), 'table went from %i to %i rows; check labels', ...
    height(t_data), height(t_annotated))
