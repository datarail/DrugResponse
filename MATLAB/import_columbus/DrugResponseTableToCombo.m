function [OutputMx, Fields, Concs] = DrugResponseTableToCombo(t_data, DrugPair, usefit)
% [OutputMx, Labels, Concs] = DrugResponseTableToCombo(t_data, DrugPair, usefit)
%   based on the columns 'DrugName' and 'DrugName2'
%
%   t_data:     table with the columns: DrugName, Conc, DrugName2, Conc2,
%                   Fields
%   usefit:     use fit for the single agent response and Bliss evaluation
%                   (default = false)
% Outputs:
%   OutputMx:    matrix for the values of Fields
%   Concs:      doses for each drugs matching Results matrix
%

checkedFields = {'CellLine' 'Time'};
assert(height(unique(t_data(:,intersect(checkedFields, varnames(t_data)))))==1)


Conc1 = unique([t_data.Conc(t_data.DrugName==DrugPair(1))
    t_data.Conc2(t_data.DrugName2==DrugPair(1))]);
Conc2 = unique([t_data.Conc(t_data.DrugName==DrugPair(2))
    t_data.Conc2(t_data.DrugName2==DrugPair(2))]);


baseFields = {'Conc' 'DrugName' 'HMSLid'};
Fields = setdiff(varnames(t_data), ...
    [checkedFields, baseFields, strcat(baseFields, '2')]);
numericField = false(length(Fields),1);
for iF = 1:length(Fields)
    numericField(iF) = isnumeric(t_data.(Fields{iF}));
end
Fields = Fields(numericField);

OutputMx = NaN(length(Conc1)+1, length(Conc2)+1, length(Fields));
if exist('usefit','var') && length(usefit)==1
    usefit = repmat(usefit, length(Fields), 1);
end

% 1st drug as a single agent
temp = t_data(t_data.DrugName==DrugPair(1) & t_data.DrugName2=='-',:);
temp = sortrows(temp , 'Conc');
[~,Cidx] = ismember(temp.Conc, Conc1);
for iF = 1:length(Fields)
    if exist('usefit','var') && usefit(iF)
        [~, ~,~,~,~,~,fit1,~,~,flag] = ICcurve_fit(temp.Conc, temp.(Fields{iF}), 'GI50');
    else flag=0;
    end
    if flag==0 % bad fit --> use the real data
        OutputMx(1+Cidx, 1, iF) = temp.(Fields{iF});
    else
        OutputMx(2:end, 1, iF) = fit1(Conc1);
    end
end

% 2nd drug as a single agent
temp = t_data(t_data.DrugName==DrugPair(2) & t_data.DrugName2=='-',:);
temp = sortrows(temp , 'Conc');
[~,Cidx] = ismember(temp.Conc, Conc2);
for iF = 1:length(Fields)
    if exist('usefit','var') && usefit(iF)
        [~, ~,~,~,~,~,fit2,~,~,flag] = ICcurve_fit(temp.Conc, temp.(varname), 'GI50');
    else flag=0;
    end
    if flag==0 % bad fit --> use the real data
        OutputMx(1, 1+Cidx, iF) = temp.(Fields{iF});
    else
        OutputMx(1, 2:end, iF) = fit2(Conc2);
    end
end

t_sub = t_data(t_data.DrugName==DrugPair(1) & t_data.DrugName2==DrugPair(2),:);
Concs = unique(t_sub.Conc2)';
for iDo = 1:length(Concs)
    temp = sortrows(t_sub(t_sub.Conc2 == Concs(iDo),:),'Conc');
    [~,Cidx] = ismember(temp.Conc, Conc1);
    for iF = 1:length(Fields)
        OutputMx(1+Cidx, 1+find(Conc2==Concs(iDo)), iF) = temp.(Fields{iF});
    end
end

Concs = {[0;Conc1] [0;Conc2]};
