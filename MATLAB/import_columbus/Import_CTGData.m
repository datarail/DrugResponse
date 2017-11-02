function t_data = Import_CTGData(folder, plateinfo, varargin)
% t_data = Import_CTGData(folder, plateinfo, varargin)
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
%addParameter(p,'T0shift', 1/4, @isscalar);
addParameter(p,'background', 0, @isvector);
addParameter(p,'TimeCourse', false, @islogical);
parse(p,varargin{:})
p = p.Results;

Files = dir(fullfile(folder, '*.xlsx'));
FileNames = {Files.name}';
Plates = unique(regexpcellsplit(FileNames,'.xlsx',1));

t_plateinfo = ImportCheckPlateInfo(plateinfo, p);
Plates = intersect(Plates, t_plateinfo.Barcode);

t_raw = [];
%ImportedFields = {};
%DefaultFields = {'Date' 'Time' 'Well'};
for iP=1:length(Plates)
    %%
    PlateFile = FileNames{regexpcell(FileNames, Plates{iP})>0};
    
    
    [~,~,temp_data] = xlsread(fullfile(folder, PlateFile));   % Read the file
    Nheader = find(cellfun(@(x)ischar(x) && any(strcmp(x,num2cell('A':'Q'))), ...
        temp_data(1:min(end,12),1)),1,'first')-1;
    assert(ismember(temp_data{Nheader+1,1},'A':'Q') && ...
        ismember(temp_data{Nheader,2},1:24)) % check the headers
    
    Cols = cellfun2(@(x)num2str(x, '%02i'), temp_data(Nheader,2:end-1));
    Rows = temp_data(Nheader+1:end,1);
    Data = cell2mat(temp_data((Nheader+1):end,2:end-1))- p.background;
    
    plate_data = table(...
        strcat(repmat(Rows,length(Cols),1), reshape(repmat(Cols,length(Rows),1),[],1)), ...
        Data(:), ...
        'VariableNames', {'Well' 'Cellcount'});
    
    plate_data = TableToCategorical(plate_data);
    

    t_raw = [t_raw; [table(repmat(Plates(iP),height(plate_data),1), ...
        'variablenames', {'Barcode'}), plate_data]];
end

t_data = AddPlateInfo_RawData(t_raw, t_plateinfo, 'Cellcount', p);
