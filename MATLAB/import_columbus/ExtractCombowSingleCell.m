function [OutputMx, OutputDist, Concs, Distbins] = ExtractCombowSingleCell(t_data, ...
    DrugPair, SingleCell, Field, Transfct, bkgrdField)
% [OutputMx, Labels, Concs] = ExtractCombowSingleCell(t_data, DrugPair, SingleCell, Variable)
%   based on the columns 'DrugName' and 'DrugName2'
%
%   t_data:     table with the columns: DrugName, Conc, DrugName2, Conc2,
%                   Fields
%   SingleCell: Corresponding SingleCell struct (indexing should match the
%                   rows in the table t_data)
%   Field:      Varible to plot (as in SinglCell); assuming the name
%                   [Field '_MeanPerWell'] for the
%   Transfct:   Function to transform the data (e.g. @log2)
%   bkgrdVar:   Variable for background subtraction (as in t_data). this is
%                   applied before the Transfct.
%
%   Data are normalized according to the control (assuming Gaussian distribution)
%
% Outputs:
%   OutputMx:   matrix for the number of cells and the mean value of
%                   Variable
%   OutputDist: SingleCell distribution for each condition and replicate
%   Concs:      doses for each drugs matching Results matrix
%   Distbins:   center of the bins for the OutputDist (same for all
%                   conditions and replicates)
%

checkedFields = {'CellLine' 'Time'};
assert(height(unique(t_data(:,intersect(checkedFields, varnames(t_data)))))==1)

if ~exist('Transfct','var')
    Transfct = @(x) x;
end

warnassert(length(SingleCell)==height(t_data), 'SingleCell is not the same length as t_data')

%%

Conc1 = unique([t_data.Conc(t_data.DrugName==DrugPair(1))
    t_data.Conc2(t_data.DrugName2==DrugPair(1))]);
Conc2 = unique([t_data.Conc(t_data.DrugName==DrugPair(2))
    t_data.Conc2(t_data.DrugName2==DrugPair(2))]);

if isvariable(t_data, 'RelGrowth')
    OutputFields = 'RelGrowth';
else
    OutputFields = 'RelCellCnt';
end
OutputFields = {OutputFields, [Field '_MeanPerWell']};
OutputMx = NaN(length(Conc1)+1, length(Conc2)+1, 2);
SCDist = cell(length(Conc1)+1, length(Conc2)+1);

% control condition
idx = find(t_data.DrugName=='-');
for iF = 1:length(OutputFields)
    OutputMx(1, 1, iF) = mean(t_data.(OutputFields{iF})(idx));
end
SCvals = [];
for i=1:length(idx)
    if ~exist('bkgrdField','var') || isempty(bkgrdField)
        SCvals = [SCvals; Transfct(SingleCell(idx(i)).(Field))];
        SCDist{1,1}{i} = SingleCell(idx(i)).(Field);
    else
        SCvals = [SCvals;
            Transfct(SingleCell(idx(i)).(Field)-t_data.(bkgrdField)(idx(i)))];
        SCDist{1,1}{i} = SingleCell(idx(i)).(Field)-t_data.(bkgrdField)(idx(i));
    end
end

Ctrl_corr = @(x) (Transfct(x)-nanmean(SCvals))/nanstd(SCvals);


% 1st drug as a single agent
for iC = 1:length(Conc1)
    idx = find(t_data.DrugName==DrugPair(1) & t_data.DrugName2=='-' & ...
        t_data.Conc==Conc1(iC));
    for iF = 1:length(OutputFields)
        OutputMx(1+iC, 1, iF) = mean(t_data.(OutputFields{iF})(idx));
    end
    for i=1:length(idx)
        if ~exist('bkgrdField','var') || isempty(bkgrdField)
            SCDist{1+iC,1}{i} = SingleCell(idx(i)).(Field);
        else
            SCDist{1+iC,1}{i} = SingleCell(idx(i)).(Field) - ...
                t_data.(bkgrdField)(idx(i));
        end
    end
end

% 2nd drug as a single agent
for iC = 1:length(Conc2)
    idx = find(t_data.DrugName==DrugPair(2) & t_data.DrugName2=='-' & ...
        t_data.Conc==Conc2(iC));
    for iF = 1:length(OutputFields)
        OutputMx(1, 1+iC, iF) = mean(t_data.(OutputFields{iF})(idx));
    end
    for i=1:length(idx)
        if ~exist('bkgrdField','var') || isempty(bkgrdField)
            SCDist{1,1+iC}{i} = SingleCell(idx(i)).(Field);
        else
            SCDist{1,1+iC}{i} = SingleCell(idx(i)).(Field) - ...
                t_data.(bkgrdField)(idx(i));
        end
    end
end

% combinations
for iC1 = 1:length(Conc1)
    for iC2 = 1:length(Conc2)
        idx = find(t_data.DrugName==DrugPair(1) & t_data.DrugName2==DrugPair(2) & ...
            t_data.Conc==Conc1(iC1) & t_data.Conc2==Conc2(iC2));

        for iF = 1:length(OutputFields)
            OutputMx(1+iC1, 1+iC2, iF) = mean(t_data.(OutputFields{iF})(idx));
        end

        for i=1:length(idx)
            if ~exist('bkgrdField','var') || isempty(bkgrdField)
                SCDist{1+iC1,1+iC2}{i} = SingleCell(idx(i)).(Field);
            else
                SCDist{1+iC1,1+iC2}{i} = SingleCell(idx(i)).(Field) - ...
                    t_data.(bkgrdField)(idx(i));
            end
        end
    end
end

Concs = {[0;Conc1] [0;Conc2]};

%% normalization
ConcIdx = ~cellfun(@isempty, SCDist);
NormSCDist = SCDist;
NormSCDist(ConcIdx) = cellfun2(@(x) cellfun2(Ctrl_corr, x), ...
    SCDist(ConcIdx));
Allvals = cellfun2(@transpose,[NormSCDist{:}]);
Allvals = [Allvals{:}];

range = [quantile(Allvals,.02) quantile(Allvals,.98)];
step = diff(range)/200;

Distbins = (step*(floor(quantile(Allvals,.02)/step)-1)):step:(step*(1+ceil(quantile(Allvals,.98)/step)));

%%
OutputDist = NormSCDist;
OutputDist(ConcIdx) = cellfun2(@(x) cellfun2(@(y) ksdensity(y,Distbins,'width',2*step)', x), ...
    NormSCDist(ConcIdx));
OutputDist(ConcIdx) = cellfun2(@(x)[x{:}],OutputDist(ConcIdx));
