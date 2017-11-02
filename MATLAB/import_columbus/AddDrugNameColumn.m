function t_out = AddDrugNameColumn(t_in, Ndrugs)

Ndrugs_in = 1+max(find(isvariable(t_in, strcat('DrugName', cellfun(@(x) {num2str(x)}, ...
    num2cell(2:9))))));

t_out = t_in;
temp_DrugName = repmat({'-'},height(t_in),1);
temp_Conc = zeros(height(t_in),1);
temp_HMSLid = temp_DrugName;
for i=(Ndrugs_in+1):Ndrugs
    t_out = [t_out table(temp_DrugName, temp_HMSLid, temp_Conc, ...
        'variableNames', strcat({'DrugName' 'HMSLid' 'Conc'}, num2str(i)))];
end
