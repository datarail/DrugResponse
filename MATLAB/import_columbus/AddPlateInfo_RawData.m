function t_data = AddPlateInfo_RawData(t_raw, t_plateinfo, NobjField, p)
% t_data = AddPlateInfo_RawData(t_raw, t_plateinfo, NobjField, p);
%
%   Incoporate the plate info to the data imported from the instrument
%
%   p are parameters of the function:
%       - p.Cellcount: function of 'NobjField' for the cell count (e.g.
%                       @(x) x(:,1)-x(:,2))
%       - p.TimeCourse
%       - p.T0date:    Input alternative to T0shift for timecourse:
%                       date and time of the treatment
%                      p.T0date is overwritten by information in the plateID file
%       - p.T0shift:   [0.25 h] Different in hours between the
%                       first plate imaged and the treatment
%

if ~isvariable(t_raw, 'Cellcount')
    t_raw = DefineCellcount(t_raw, NobjField, p);
end

plate_barcodes = unique(t_raw.Barcode);

warnassert(length(plate_barcodes)>=height(t_plateinfo), ...
    'Found %i plates, but expected %i plates; missing %s', ...
    length(plate_barcodes), height(t_plateinfo), ...
    strjoin(setdiff(t_plateinfo.Barcode,plate_barcodes)',' ') );

otherVariables = setdiff(varnames(t_plateinfo), {'Time' 'CellLine' 'Barcode' ...
    'TreatmentFile' 'DesignNumber' 'ExpNumber'}, 'stable');
if ismember('Untrt', otherVariables)
    if ~all(strcmp(t_plateinfo.TreatmentFile(t_plateinfo.Untrt>0),'-')) || ...
            ~all(t_plateinfo.Untrt(strcmp(t_plateinfo.TreatmentFile,'-'))>0)
        warnprintf('Discrepency between TreatmentFile==''-'' and Untrt; Overwriting TreatmentFile')
        t_plateinfo.TreatmentFile(t_plateinfo.Untrt>0) = {'-'};
    end
    t_plateinfo.Untrt = [];
    otherVariables = setdiff(otherVariables, 'Untrt', 'stable');
end

CellLine = cell(height(t_raw),1);
Barcode = cell(height(t_raw),1);
TreatmentFile = cell(height(t_raw),1);
DesignNumber = zeros(height(t_raw),1);
Time = zeros(height(t_raw),1);
Untrt = false(height(t_raw),1);


cnt = 0;
% join the barcode table to t_raw with a few controls
for iBC = 1:height(t_plateinfo)

    idx = find(strcmp(t_raw.Barcode, t_plateinfo.Barcode{iBC}));
    if isempty(idx)
        idx = find(strfindcell(t_raw.Barcode, t_plateinfo.Barcode{iBC}));
        if isempty(idx)
            warnprintf('No result found for plate %s', t_plateinfo.Barcode{iBC})
        else
            warnprintf('Partial barcode match for plate %s', t_plateinfo.Barcode{iBC})
        end
    end

    CellLine(idx) = t_plateinfo.CellLine(iBC);
    Barcode(idx) = t_plateinfo.Barcode(iBC);
    TreatmentFile(idx) = t_plateinfo.TreatmentFile(iBC);
    [~,~,ext] = fileparts(t_plateinfo.TreatmentFile{iBC});

    if (strcmp(ext,'.mat') || strcmp(ext,'.hpdd'))
        assert(isvariable(t_plateinfo,'DesignNumber'), ...
            'Treatment files .mat or .hpdd needs a DesignNumber column')
        DesignNumber(idx) = t_plateinfo.DesignNumber(iBC);
    elseif ~isempty(ext) % assume a tab-separated file
        if isvariable(t_plateinfo,'DesignNumber') && ...
                (t_plateinfo.DesignNumber(iBC)~=1 || ...
                ~isnan(t_plateinfo.DesignNumber(iBC)))
            warnprintf('for tab-separated file, DesignNumber is forced to be 1')
        end
        DesignNumber(idx) = 1;
    end
    Untrt(idx) = strcmp(t_plateinfo.TreatmentFile(iBC),'-');
    if isfield(p, 'TimeCourse') && p.TimeCourse
        % rounded to the 1/100 of hour, Time = (scan Time) - T0date
        % if no T0date given, 1st value is p.T0shift (default ~15 minutes)
        if isvariable(t_plateinfo, 'T0date')
            Time(idx) = .01*round(100*( 24*(t_raw.Date(idx)-datenum(t_plateinfo.T0date(iBC)) )));
        elseif ~isempty(p.T0date)
            Time(idx) = .01*round(100*( 24*(t_raw.Date(idx)-datenum(p.T0date) )));
        else
            Time(idx) = .01*round(100*( 24*(t_raw.Date(idx)-min(t_raw.Date)) + p.T0shift ));
        end
    else
        temp = t_plateinfo.Time(iBC);
        if iscellstr(temp); temp = cellstr2mat(temp); end
        Time(idx) = temp;
    end

    % parse the additional plate information from the barcode file
    for i = 1:length(otherVariables)
        eval([otherVariables{i} '(idx) = t_plateinfo.' otherVariables{i} '(iBC);'])
    end

    cnt = cnt+1;
end

% format properly the additional plate information from the plateinfo file
for i = 1:length(otherVariables)
    eval([otherVariables{i} ' = ToColumn(' otherVariables{i} ');'])
    eval(['if length(' otherVariables{i} ')<height(t_raw), ' ...
        otherVariables{i} '(height(t_raw)+1) = ' otherVariables{i} '(1);' ...
        otherVariables{i} '(height(t_raw)+1) = []; end'])
end

% broadcast the properties
if cnt<length(plate_barcodes)
    warning('some entries in the result file are unused!')
    Usedidx = ~cell2mat(cellfun2(@isempty,CellLine));
    TreatmentFile = TreatmentFile(Usedidx);
    CellLine = CellLine(Usedidx);
    DesignNumber = DesignNumber(Usedidx);
    Time = Time(Usedidx);
    Barcode = Barcode(Usedidx);
    for i = 1:length(otherVariables)
        eval([otherVariables{i} ' = ' otherVariables{i} '(Usedidx);'])
    end
else
    Usedidx = 1:height(t_raw);
end

% Evaluate the untreated plates
Untrt = cellfun(@(x) strcmp(x,'-') || isempty(x), TreatmentFile) ;
assert(~any(DesignNumber==0 & ~Untrt), 'Some wells are not ''Untrt'' and don''t have a DesignNumber')
warnassert(all(Time(~Untrt)>0), 'Some treated/perturbed wells have Time=0')

% compile the final table
t_data = [table(Barcode, CellLine, TreatmentFile, DesignNumber, Untrt, Time) ...
    t_raw(Usedidx, intersect([{'Well' 'Date' 'Cellcount' 'Deadcount'} ToRow(NobjField)], varnames(t_raw), 'stable'))];
if ~isempty(otherVariables)
    for i = 1:length(otherVariables)
        if isvariable(t_data, otherVariables{i})
            fprintf(['\tReplacing variable(s): ' otherVariables{i} ' by value in Plate info file\n'])
            eval(['t_data = [t_data(:,setdiff(varnames(t_data),''' otherVariables{i} ''',''stable''))'...
                ' table(' otherVariables{i} ')];'])
        else
            fprintf(['\tAdding variable(s): ' otherVariables{i} ' from Plate info file\n'])
            eval(['t_data = [t_data'...
                ' table(' otherVariables{i} ')];'])
        end
    end
end

t_data = TableToCategorical(t_data);

fprintf('\n')
