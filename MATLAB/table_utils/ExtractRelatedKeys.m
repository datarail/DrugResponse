function keySets = ExtractRelatedKeys(variableNames, keys)
% keySets = ExtractRelatedKeys(variableNames, keys)

keySets = cell(1,length(keys));
for i=1:length(keys)
    keySets{i} = variableNames(regexpcell(variableNames, ['^' keys{i}])==1);
end
