
function [mRNA_UMI, mRNA_labels, mRNA_total] = ImportAnnotate_DGEdata(filename, plateinfo, Barcode, folder)
% [mRNA_values, mRNA_labels] = ImportAnnotate_DGEdata(filename, plateinfo, Barcode, folder)

if iscellstr(filename)
    assert(length(filename)==length(Barcode))
    
    for i=1:length(filename)
        if nargout==3
            [mRNA_UMI{i}, mRNA_labels{i}, mRNA_total{i}] = ...
                ImportAnnotate_DGEdata(filename{i}, plateinfo, Barcode{i}, folder);
        else
            [mRNA_UMI{i}, mRNA_labels{i}] = ...
                ImportAnnotate_DGEdata(filename{i}, plateinfo, Barcode{i}, folder);
        end
    end
    n_cond = sum(cellfun(@(x) size(x,2), mRNA_UMI));
    for i=2:length(filename)
        assert(all(mRNA_labels{1}{1}.gene==mRNA_labels{i}{1}.gene))
        mRNA_UMI{1} = [mRNA_UMI{1} mRNA_UMI{i}];
        mRNA_labels{1}{2} = [mRNA_labels{1}{2}; mRNA_labels{i}{2}];
        if exist('mRNA_total', 'var')
            mRNA_total{1} = [mRNA_total{1} mRNA_total{i}];
        end
    end
    assert(n_cond==size(mRNA_UMI{1},2));
    assert(n_cond==height(mRNA_labels{1}{2}));
    mRNA_UMI = mRNA_UMI{1};
    mRNA_labels = mRNA_labels{1};
    if exist('mRNA_total', 'var')
        mRNA_total = mRNA_total{1};
    end
    return
end

%%
p = struct('TimeCourse', false, 'Cellcount', 'none');
t_plateinfo = ImportCheckPlateInfo(plateinfo, p);


RNAfilename = [filename '.unq.refseq.umi.dat'];
rawRNA = importdata(RNAfilename);


Wells = CorrectWellsTo0xY(regexpcelltokens(rawRNA.textdata(1,2:end), '_([A-Q][0-9]*)'))';
t_labels = table(repmat({Barcode}, length(Wells),1), Wells, ...
    'variablenames', {'Barcode' 'Well'});

t_labels = AddPlateInfo_RawData(t_labels, t_plateinfo, '', p);
t_labels = Annotate_CellCountData(t_labels, folder);

INFOfilename = [filename '.unq.well_summary.dat'];
rawINFO = importdata(INFOfilename);
Wells = CorrectWellsTo0xY(regexpcelltokens(rawRNA.textdata(1,2:end), '_([A-Q][0-9]*)'))';
t_info = table(categorical(Wells), ...
    rawINFO.data(find(strcmp(rawINFO.textdata(:,1),'Refseq_Total'))-1,:)', ...
    rawINFO.data(find(strcmp(rawINFO.textdata(:,1),'Refseq_UMI'))-1,:)', ...
    'variablenames', {'Well' 'Refseq_Total' 'Refseq_UMI'});

t_labels = leftjoin(t_labels, t_info);

mRNA_names = table(categorical(rawRNA.textdata(2:end,1)),'VariableNames',{'gene'});
mRNA_UMI = rawRNA.data;

mRNA_labels = {mRNA_names t_labels};

if nargout==3
    RNAfilename = [filename '.unq.refseq.total.dat'];
    rawRNA = importdata(RNAfilename);
    assert(all(mRNA_names.gene==rawRNA.textdata(2:end,1)))
    assert(all(t_info.Well==...
        CorrectWellsTo0xY(regexpcelltokens(rawRNA.textdata(1,2:end), '_([A-Q][0-9]*)'))'));
    mRNA_total = rawRNA.data;
end
