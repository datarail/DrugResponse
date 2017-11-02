function t_data = Import_PlatesCellCountData2(filename, plateinfo, varargin)
% t_data = Import_PlatesCellCountData2(filename, plateinfo, varargin)
%
%   merge the information in the barcode (either file or table) with the
%   data from the output of Columbus or from other tsv file with proper
%   columns.
%
%   filename is a tsv file which should contain at least the headers
%   (default output from Columbus):
%           - Result (containing the barcodes and the timestamp)
%               -> or two columns Barcode and Date
%           - Well
%           - NumberOfAnalyzedFields
%           - Nuclei - Number Of Objects; can be changed to another field with
%               option 'NobjField', e.g.
%               Import...(..., 'NobjField', 'Nuclei Selected - Number of Objects Nuclei Selected')
%           note: this allows to get multiple fields from the raw datafile;
%               the first field will be converted to Cellcount, other field
%               names will not be changed (except being made compatible).
%               For example: Import...(..., 'NobjField', ...
%                   {'Nuclei Selected - Number of Objects', ...
%                    'Nuclei - Number Of Objects'} )
%
%   Alternatively, filename can be a Nx[1-2] cell array with a filename and
%   corresponding barcode. Files will be read in order and if it is a Nx2
%   array, barcode will be added while reading. Each file should have
%   fields Well, NumberOfAnalyzedFields, Nuclei_NumberOfObjects, and a
%   field Result/Barcode if a Nx1 array is provided.
%
%   plateinfo should be a table (or the name of a tsv file) with headers:
%           - Barcode
%           - CellLine
%           - TreatmentFile (refer to a table or .hpdd/.mat with with the
%               plate design; use  -  if untreated, 1 if tab-separated file)
%           - DesignNumber (replicate or treatment design, mandatory for
%               .mat or .hpdd treatment files; set to 1 for other files)
%           - Replicate (optional - to differentiate to plates with the same
%               treatment file, and design number if .mat/.hpdd)
%           - Time (in hours);
%                   note: Time is ignored if the option 'TimeCourse' is
%                   true; in that case, the timestamp in the Columbus file
%                   will be used to tag the time of the sample.
%           - any Additional field with relevant plate properties to pass to the
%               data (collagen, ...)
%           note: 'ExpNumber' is ignored as it serves as an internal
%               control for the D300
%
%   varargin:   - 'NobjField'   [Nuclei - Number Of Objects]
%               - 'Cellcount'   [Nuclei - Number Of Objects]
%               - 'TimeCourse'  [false]
%               - 'T0shift'     [0.25 h] Different in hours between the
%                                   first plate imaged and the treatment
%               - 'T0date'      Input alternative to T0shift for timecourse:
%                                   date and time of the treatment
%
%   t_data is a table with each well annotated accorting to the barcode.
%   The column 'Untrt' is evaluated and the data are corrected for the
%   number of fields (stored in column Cellcount).
%       the column Cellcount can be a function of the input columns (by
%       default Cellcount is 'Nuclei_NumberOfObjects' or the first column
%       in the field 'NobjField').
%
%   Example:
%   --------
%
%   t_data = Import_PlatesCellCountData('Results_20150320.txt'], ...
%     'PlateIDs_20150320.txt', 'NobjField', ...
%     {'Nuclei-Hoechst - Number of Objects' 'Nuclei-LDRpos - Number of Objects'},...
%     'Cellcount', @(x) x(:,1)-x(:,2), 'Deadcount', @(x) x(:,2));
%

fprintf('Import Cell count data from Columbus files:\n');

p = inputParser;
addParameter(p,'NobjField',{'Nuclei_NumberOfObjects'},@(x) ischar(x) || iscellstr(x));
addParameter(p,'Cellcount', [], @(x) isa(x,'function_handle'));
addParameter(p,'Deadcount', [], @(x) isa(x,'function_handle'));
addParameter(p,'TimeCourse', false, @islogical);
addParameter(p,'RemoveNobjField', true, @islogical);
addParameter(p,'T0shift', 1/4, @isscalar);
addParameter(p,'T0date', [], @(x) ~isnan(datenum(x)));
parse(p,varargin{:})
p = p.Results;

NobjField = p.NobjField;
if ischar(NobjField), NobjField = {NobjField}; end
try
    NobjField = matlab.internal.tableUtils.makeValidName(NobjField,'silent');
catch
    NobjField = ReducName( ReplaceName(NobjField,['-/\:?!,.' 250:1e3], '_'), ' ');
end

%% check for proper inputs
if ischar(filename)
    assert(exist(filename,'file')>0, 'File %s is missing!', filename)
else
    assert(iscellstr(filename), 'Wrong type of input')
    for i=1:size(filename,1)
        assert(exist(filename{i,1},'file')>0, 'File %s is missing!', filename{i,1})
    end
end


%% load the cell count data

if ischar(filename)
    % case of a single file
    fprintf('\nImporting from: %s\n', filename)
    t_raw = tsv2table(filename);
else
    % case of a cell array (multiple files)
    t_raw = table;
    fprintf('\nImporting from: files\n')
    for i=1:size(filename,1)
        temp = tsv2table(filename{i,1});
        fprintf('\t%s\n', filename{i,1});

        if size(filename,2)==2
            % replace/ add the barcode
            if any(isvariable(temp, {'Result' 'Barcode'}))
                warnprintf(['Overwriting the Result/Barcode found in ' ...
                    'the individual file with input filename(:,2)'])
                if isvariable(temp, 'Barcode')
                    temp(:,'Barcode') =[];
                else
                    temp(:,'Result') =[];
                end
            end
            temp = [table(repmat(filename(i,2), height(temp),1), ...
                'VariableName', {'Barcode'}) temp];
        end

        t_raw = [t_raw; temp];

    end
end


% specific case of the output of Columbus: Results split into Barcode and
% date
if isvariable(t_raw, 'Result')
    t_raw = splitBarcodeDate(t_raw);
end

% correction of some variable names
CorrectionVarNames = {
    'WellName' 'Well'
    'PlateName' 'Barcode'
    'MeasurementDate' 'Date'};
for i=1:size(CorrectionVarNames,1)
    if isvariable(t_raw, CorrectionVarNames{i,1})
        t_raw.Properties.VariableNames{CorrectionVarNames{i,1}} = ...
            CorrectionVarNames{i,2};
    end
end

% report error if missing field for cell count
if ~ismember(NobjField{1}, varnames(t_raw))
    warnprintf('Missing the field %s used for cell count', NobjField{1})
    disp('Available fields:')
    disp(varnames(t_raw)')
    error('Use optional input ''NobjField'' for specifying fields')
end
if length(NobjField)>1 && ~strcmp(NobjField{2},'all')
    if ~all(ismember(NobjField, varnames(t_raw)))
        warnprintf('Not all fields in ''NobjField'' are present in the file')
        disp('Available fields:')
        disp(varnames(t_raw)')
        error('Missing fields: %s', strjoin(NobjField(~ismember(NobjField, varnames(t_raw)))))
    end
elseif length(NobjField)==2 && strcmp(NobjField{2},'all')
    NobjField = [NobjField(1) setdiff(varnames(t_raw), [NobjField(1) ...
        {'Well' 'Barcode' 'Result' 'Date' 'Row' 'Column' ...
        'NumberOfAnalyzedFields' 'URL'}])];
end

% check the number of fields
if length(unique(t_raw.NumberOfAnalyzedFields))>1
    Nref = median(t_raw.NumberOfAnalyzedFields);
    warnprintf('%i wells with missing fields', sum(t_raw.NumberOfAnalyzedFields<Nref))
    FieldCorrected = NobjField(strfindcell(NobjField,'Number')>0);
    warnprintf('Correcting for field number:\n\t - %s', strjoin(FieldCorrected,'\n\t - '))
    for i = find(t_raw.NumberOfAnalyzedFields~=Nref)'
        for j=1:length(FieldCorrected)
            t_raw.(FieldCorrected{j})(i) = t_raw.(FieldCorrected{j})(i)*...
                (Nref/t_raw.NumberOfAnalyzedFields(i));
        end
    end
end

t_raw.Well = CorrectWellsTo0xY(t_raw.Well);

%%

try
    t_plateinfo = ImportCheckPlateInfo(plateinfo, p);
catch err
    warnprintf('Error with the plate info file (from function ImportCheckPlateInfo):')
    warning(err.message)
    warnprintf(' --> Annotation discarded; output is raw data')
    t_data = TableToCategorical(DefineCellcount(t_raw, NobjField, p));
    return
end
t_data = AddPlateInfo_RawData(t_raw, t_plateinfo, NobjField, p);

if p.RemoveNobjField
    fprintf('Removing field: %s\n', NobjField{:});
    for iN = 1:length(NobjField)
        t_data.(NobjField{iN}) = [];
    end
end

fprintf('\n')
end


function t_raw = splitBarcodeDate(t_raw)
Code_date = regexp(t_raw.Result,' > ','split');
Code_date = vertcat(Code_date{:});
Code_date(:,2) = regexpcellsplit(Code_date(:,2), '(',1);
Code_date(:,2) = cellfun2(@datenum, Code_date(:,2));
t_raw = [cell2table(Code_date,'VariableName',{'Barcode' 'Date'}) t_raw(:,2:end)];
end
