function t_data = Import_IncucyteCellCountData(folder, plateinfo, varargin)
% t_data = Import_IncucyteCellCountData(folder, plateinfo, varargin)
%
%   Merge the information in the barcode (either file or table) with the
%   data from the output of the Incucyte or from other tsv file with proper
%   columns.
%   Files should be all in the 'folder' and the first part of the filename
%   (without space) is the barcode of the plate reported in the
%   'plateinfo', the second part of the filename is the 'field'
%
%   varargin:   - 'NobjField'   [total cells]
%               - 'Cellcount'   [total cells]
%               - 'T0shift'     [0.25 h] Different in hours between the
%                                   first plate imaged and the treatment
%               - 'T0date'      Input alternative to T0shift for timecourse:
%                                   date and time of the treatment
%
%   t_data is a table with each well annotated accorting to the barcode.
%   The column 'Untrt' is evaluated and the data are corrected for the
%   number of fields (stored in column Cellcount).
%       the column Cellcount can be a function of the input columns (by
%       default Cellcount is 'total cells' or the first column
%       in the field 'NobjField').


p = inputParser;
addParameter(p,'NobjField',{'total cells'},@(x) ischar(x) || iscellstr(x));
addParameter(p,'Cellcount', [], @(x) isa(x,'function_handle'));
addParameter(p,'T0shift', 1/4, @isscalar);
addParameter(p,'T0date', [], @isvector);
parse(p,varargin{:})
p = p.Results;
p.TimeCourse = true;

Files = dir(fullfile(folder, '*.txt'));
FileNames = {Files.name}';
Plates = unique(regexpcellsplit(FileNames,' ',1));

NobjField = p.NobjField;
if ischar(NobjField), NobjField = {NobjField}; end
try
    NobjField = matlab.internal.tableUtils.makeValidName(NobjField,'silent');
catch
    NobjField = ReducName( ReplaceName(NobjField,['-/\:?!,.' 250:1e3], '_'), ' ');
end

t_plateinfo = ImportCheckPlateInfo(plateinfo, p);

t_raw = [];
ImportedFields = {};
DefaultFields = {'Date' 'Time' 'Well'};
for iP=1:length(Plates)
    %%
    PlateFiles = FileNames(regexpcell(FileNames, Plates{iP})>0);
    plate_data = cell(length(PlateFiles),1);
    for iF = 1:length(PlateFiles)
        %%
        field = regexp(PlateFiles{iF}, '\ ([\w\ \d]*).txt', 'tokens');
        field = field{1};
        try
            fieldname = matlab.internal.tableUtils.makeValidName(field, 'silent');
        catch
            fieldname = ReducName( ReplaceName(fieldname,['-/\:?!,.' 250:1e3], '_'), ' ');
        end

        temp_data = tsv2cell(fullfile(folder, PlateFiles{iF}));   % Read the file
        Nheader = find(strcmp(temp_data(:,1), 'Date Time')); % get the position of data
        assert(length(Nheader)==1)
        assert(strcmp(temp_data(Nheader,2), 'Elapsed')) % check the headers
        assert(ismember(temp_data{Nheader,3}(1),'A':'Q') & ...
            ismember(temp_data{Nheader,3}(2),'0':'9')) % check the headers

        Date = cellfun(@datenum, temp_data((Nheader+1):end,1));
        Time = cellfun(@str2double, temp_data((Nheader+1):end,2));
        Wells = CorrectWellsTo0xY(temp_data(Nheader,3:end));
        Data = cellfun(@str2double,temp_data((Nheader+1):end,3:end));

        plate_data{iF} = table(repmat(Date, length(Wells),1), repmat(Time, length(Wells),1), ...
            reshape(repmat(Wells,length(Date),1),[],1), Data(:), ...
            'VariableNames', [DefaultFields fieldname]);

        plate_data{iF} = TableToCategorical(plate_data{iF});

    end

    if length(PlateFiles)>1
        assert(length(unique(cellfun(@height,plate_data)))==1)

        for iF = 2:length(PlateFiles)
            plate_data{1} = innerjoin(plate_data{1}, plate_data{iF}, 'keys', ...
                DefaultFields);
        end
        assert(height(plate_data{1}) == height(plate_data{2}))
        assert(width(plate_data{1}) == (length(DefaultFields)+length(PlateFiles)))
    end

    if iP==1
        ImportedFields = setdiff(varnames(plate_data{1}),DefaultFields);
    else
        if length(varnames(plate_data{1}))~=(length(ImportedFields)+length(DefaultFields)) || ...
                ~all(ismember(ImportedFields, varnames(plate_data{1})))
            disp(['Imported fields: ' strjoin(ImportedFields,'; ')])
            disp(['fields for ', Plates{iP} ': ' ...
                strjoin(plate_data{1}.Properties.VariableNames((length(DefaultFields)+1):end),'; ')])
            error('Mismatch between import fields and fields for %s', Plates{iP})
        end
    end

    t_raw = [t_raw; [table(repmat(Plates(iP),height(plate_data{1}),1), ...
        'variablenames', {'Barcode'}), plate_data{1}]];
end

assert(all(ismember(NobjField, ImportedFields)), 'Some specified fields were not found')
if any(~ismember(ImportedFields, NobjField))
    fprintf('Ignored field: %s\n', ImportedFields{~ismember(ImportedFields, NobjField)})
end

t_data = AddPlateInfo_RawData(t_raw, t_plateinfo, NobjField, p);
