
function [SingleCelldata, t_SingleCelldata] = Import_CyCIFData(t_conditions, ...
    folder, subfolder, IF_key_file, process_key_file)
% [SingleCelldata, t_SingleCelldata] = Import_CyCIFData(t_conditions, ...
%         folder, subfolder, IF_key_file, process_key_file)

defaultFields = {
    'Nuclei - Intensity Nucleus Alexa %s Mean'
    'Nuclei - N/C %s'
    'Nuclei - Intensity Nucleus HOECHST 33342 Mean'
    'Nuclei - Nucleus Area [µm²]'
    'Nuclei - DNA content'
    'X'
    'Y'
    'Field'};

IFkey = tsv2table([folder IF_key_file]);

cycles = varnames(IFkey);
t_Abs = IFkey(:,['Barcode' cycles(strfindcell(cycles, 'cycle')==1)]);

t_plate_keys = tsv2table([folder process_key_file]);
cycle_data = setdiff(varnames(t_plate_keys),'Barcode');


% check that the file exists
for iP = 1:height(t_plate_keys)
    for iR = 1:length(cycle_data)
        f = dir([folder subfolder t_plate_keys.Barcode{iP} '/*' ...
            t_plate_keys.(cycle_data{iR}){iP} '*txt']);
        f = {f.name}';
        f = f(regexpcell(f, 'Whole Image')==0);
        assert(length(f)>=sum(t_conditions.Barcode==t_plate_keys.Barcode{iP}), ...
            'files missing for %s, %s (tag = |%s|)', t_plate_keys.Barcode{iP}, ...
            cycle_data{iR}, t_plate_keys.(cycle_data{iR}){iP})
    end
end

for iP = 1:height(t_plate_keys)
    subt_conditions = t_conditions(t_conditions.Barcode==t_plate_keys.Barcode{iP},:);

    for iR = 1:length(cycle_data)
        fprintf('\n---------------------\n%s: %s\n-------------------\n',...
            t_plate_keys.Barcode{iP}, cycle_data{iR});

        Abs = t_Abs(strcmp(t_Abs.Barcode,t_plate_keys.Barcode(iP)), ...
            strfindcell(varnames(t_Abs), cycle_data{iR})==1);
        channel = cellstr2mat(regexpcelltokens(varnames(Abs), 'cycle[0-9]_([0-9]*)'));
        channel_idx = ~strcmp([Abs{1,:}], '-') & ~isempty([Abs{1,:}]);
        channel = channel(channel_idx);
        Abs = [Abs{1,channel_idx}];


        SelectedFields = {};
        FieldName = {};
        for iF = 1:length(defaultFields)
            if ~isempty(regexp(defaultFields{iF}, '%s', 'once'))
                for iC = 1:length(channel)
                    SelectedFields{end+1,1} = sprintf(defaultFields{iF}, num2str(channel(iC)));
                    FieldName{end+1,1} = sprintf([defaultFields{iF} '_%s'], ['_' Abs{iC} '_'], cycle_data{iR});
                end
            else
                SelectedFields{end+1,1} = defaultFields{iF};
                FieldName{end+1,1} = sprintf([defaultFields{iF} '_%s'], cycle_data{iR});
            end
        end
        FieldName = matlab.internal.tableUtils.makeValidName(FieldName, 'silent');
        SelectedFields = matlab.internal.tableUtils.makeValidName(SelectedFields, 'silent');

        temp_SingleCelldata = Import_SingleCell(subt_conditions, ...
            [folder subfolder], SelectedFields, false, t_plate_keys.(cycle_data{iR}){iP});

        for iF = 1:length(temp_SingleCelldata)
            for iC = 1:length(channel)
                temp_SingleCelldata(iF).Background.Properties.VariableNames = ...
                    matlab.internal.tableUtils.makeValidName( ...
                    regexprep(varnames(temp_SingleCelldata(iF).Background), ...
                    num2str(channel(iC)), ['_' Abs{iC} '_' cycle_data{iR}]), 'silent');
            end
        end

        for iF = 1:length(FieldName)
            [temp_SingleCelldata.(FieldName{iF})] = temp_SingleCelldata.(SelectedFields{iF});
            temp_SingleCelldata = rmfield(temp_SingleCelldata, SelectedFields{iF});
        end

        if iR==1
            Allcycles = temp_SingleCelldata;
        else
            fields = setdiff(fieldnames(temp_SingleCelldata),'Background');
            Allfields = fieldnames(Allcycles);
            for iF = 1:length(temp_SingleCelldata)
                Allcycles(iF).Background = leftjoin(Allcycles(iF).Background, ...
                    temp_SingleCelldata(iF).Background, 'Keys', 'Field');
            end
            for iF = 1:length(fields)
                if ~ismember(fields{iF}, Allfields)
                    [Allcycles.(fields{iF})] = temp_SingleCelldata.(fields{iF});
                end
            end
        end
    end

    if iP==1
        SingleCelldata = Allcycles;
        t_SingleCelldata = subt_conditions;
    else
        SingleCelldata = [SingleCelldata; Allcycles];
        t_SingleCelldata = [t_SingleCelldata; subt_conditions];
    end
end

fprintf('\n---------------------\n done\n\n');
