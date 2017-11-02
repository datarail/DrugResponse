function [SingleCell, t_out] = Import_SingleCell(t_processed, folder, fields, timecourse, filefilter)
% [SingleCell, t_out] = Import_SingleCell(t_processed, folder, fields, timecourse, filefilter)
%   t_processed     table outputed from Merge_CellCountData; need the
%                       columns 'Barcode' 'Well' and will condiser the data
%                       columns that ends with '_MeanPerWell' (otherwise,
%                       specify 'fields' to select columns in the single
%                       cell data).
%                       If column 'Date' exists, it will assume a timecourse
%                       and look for corresponding file name (unless
%                       variable timecourse is false)
%
%   folder          folder in which the single cell data are stored
%
%   fields          fields to kepp in the single cell data
%
%   timecourse      time course experiment, multiple time points for the
%                       same well (optional, default=false)
%
%   filefilter      regular experssion for filtering files in case of
%                       multiple output files for the same plate/well
%
%   SingleCell      structure with the data fields for each plate/well
%
%   t_out           table with properties for the SingleCell (equal to
%                       t_processed unless multiple files/times were found)
%

%%
if ~exist('timecourse', 'var') || isempty(timecourse)
    timecourse = isvariable(t_processed,'Date');
end

if ~exist('filefilter', 'var')
    filefilter = '';
end

if isempty(fields)
    f = ls([folder '/' char(t_processed.Barcode(1)) '/*' filefilter '*txt']);
    ftemp = fopen([folder '/' char(t_processed.Barcode(1)) '/' f(1,:)],'r');
    fields = regexp(fgetl(ftemp),'\t','split');
    warnprintf('Available fields are:\n\t - %s \n', strjoin(fields,'\n\t - '));
    error('no field found')
else
    try
        fields = matlab.internal.tableUtils.makeValidName(fields,'silent');
    catch
        fields = ReducName( ReplaceName(fields,['-/\:?!,.' 250:1e3], '_'), ' ');
    end
    fprintf('\nFields are:\n - %s\n', strjoin(fields,'\n - '));
end

SingleCell = struct;
for i=1:length(fields)
    SingleCell.(fields{i}) = [];
end
SingleCell.Barcode = ''; SingleCell.Well = ''; SingleCell.Date = '';
SingleCell.Background = [];
SingleCell = repmat(SingleCell, height(t_processed),1);
%%

assert(exist(folder,'dir')>0, ['Folder: ' folder ' missing'])

if isvariable(t_processed,'Date')
    t_out = t_processed;
else
    t_out = [t_processed table(repmat(0, height(t_processed), 1), ...
        'variablenames', {'Date'})];
end
output_cnt = 0;

fprintf('\nReading all conditions (%i total):\n', height(t_processed));
for iPW = 1:height(t_processed)
    loop_waitbar(iPW, height(t_processed))

    folderlist = dir(folder);
    allfolder = setdiff({folderlist([folderlist.isdir]==1).name}, {'.' '..'});

    Plate = char(t_processed.Barcode(iPW));
    Well = char(t_processed.Well(iPW));

    if Well(2)=='0', Well = Well([1 3]); end


    subfolder = [folder filesep allfolder{regexpcell(allfolder, Plate)>0}];

    if exist('filefilter','var')
        files = dir([subfolder filesep '*' filefilter '*.t*']); % for tsv or txt files
        if isempty(files)
            dirs = dir([subfolder filesep '*']);
            dirs = setdiff({dirs([dirs.isdir]).name},{'.', '..'});
            subfolder = [subfolder filesep dirs{1}];
            files = dir([subfolder filesep '*' filefilter '*.t*']);
        end
    else
        files = dir([subfolder filesep '*.t*']); % for tsv or txt files
    end
    files = {files([files.isdir]==0).name};
    files = files(regexpcell(files, sprintf('result.%s\\[', Well))>0);

    if isvariable(t_processed,'Date') && timecourse
        files = files(regexpcell(files, ...
            datestr(t_processed.Date(iPW),'yyyy-mm-ddTHHMMSS'))==1);
    end


    if isempty(files)
        warnprintf(['Missing file for: ' Plate ' ' Well]);
        continue
    end

    files_WI = files(regexpcell(files, 'Whole Image')>0);
    files = files(regexpcell(files, 'Whole Image')==0);
    if ~timecourse || isvariable(t_processed,'Date')
        warnassert(length(files)==1, ['Too many files: ' strjoin(files,' || ')])
        files = files(1);
    end

    for it = 1:length(files)
        try
            date = regexp(files{it},'^([0-9\-]*)T','tokens');
            time = regexp(files{it}, '^[0-9\-]*T([0-9]*)\-', 'tokens');
            Date = datenum([date{1}{1} '-' time{1}{1}], 'yyyy-mm-dd-HHMMSS');
        catch
            Date = '';
        end

        t_ss = tsv2table([subfolder filesep files{it}]);

        output_cnt = output_cnt+1;

        if isvariable(t_processed,'Date')
            t_out(output_cnt,:) = t_processed(iPW,:);
        else
            t_out(output_cnt,:) = [t_processed(iPW,:) table(Date)];
        end
        if ~isempty(files_WI)
            t_bkgrd = tsv2table([subfolder filesep files_WI{it}]);
            bck_fields = varnames(t_bkgrd);
            SingleCell(output_cnt).Background = t_bkgrd(:,['Field' ...
                bck_fields(regexpcell(bck_fields, '[0-9]*Mean')>0)]);
        end

        SingleCell(output_cnt).Barcode = Plate;
        SingleCell(output_cnt).Well = Well;
        SingleCell(output_cnt).Date = Date;
        for i=1:length(fields)
            if ~isvariable(t_ss, fields{i})
                ftemp = fopen([subfolder filesep files{it}],'r');
                allfields = regexp(fgetl(ftemp),'\t','split');
                warnprintf('Available fields are:\n\t - %s \n', strjoin(allfields,'\n\t - '));
                error('Field %s not found in file %s', fields{i}, files{it})
            end
            SingleCell(output_cnt).(fields{i}) = t_ss.(fields{i});
        end
    end

end

fprintf('\n');
