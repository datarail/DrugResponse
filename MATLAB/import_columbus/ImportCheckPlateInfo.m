function t_plateinfo = ImportCheckPlateInfo(plateinfo, p)

if ischar(plateinfo)
    assert(exist(plateinfo,'file')>0, 'Barcode file %s is missing!', plateinfo)
    t_plateinfo = tsv2table(plateinfo);
elseif istable(plateinfo)
    t_plateinfo = plateinfo;
else
    error('Wrong argument for plateinfo')
end

CheckPlateInfo(t_plateinfo, p.TimeCourse)
if isvariable(t_plateinfo, 'DesignNumber') && iscell(t_plateinfo.DesignNumber)
    t_plateinfo.DesignNumber = cellstr2mat(t_plateinfo.DesignNumber);
end

end


function CheckPlateInfo(t_plateinfo, IsTimeCourse)
Infovars = {'Barcode' 'Time' 'CellLine' 'TreatmentFile'};
if IsTimeCourse, Infovars = setdiff(Infovars, 'Time'); end
for i=1:length(Infovars)
    assert(isvariable(t_plateinfo, Infovars{i}), 'Missing columns %s in the plate info',...
        Infovars{i})
end
for i=1:height(t_plateinfo)
    tf = t_plateinfo.TreatmentFile{i};
    [~,n,e] = fileparts(tf);
    assert(ismember(e, {'.mat' '.hpdd' '.tsv' '.txt'}) | strcmp(n,'-'), ...
        'TreatmentFile should be a .mat, .hpdd, .txt, .tsv file of ''-'' for untreated')
end
end
