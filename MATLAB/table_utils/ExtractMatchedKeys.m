function Dimensions = ExtractMatchedKeys(t_data, keys)
% Dimensions = ExtractMatchedKeys(t_data, keys)

keySets = cell(1,0);
Dimensions = keySets;
dataVars = varnames(t_data);
assert(all(ismember(keys, dataVars)), ['%s are not variables of input table' ...
    '; variables are: \n\t%s'], ...
    strjoin(keys(~ismember(keys, dataVars)),', '), ...
    strjoin(dataVars, '\n\t'))

for i=1:length(keys)
    % check that the key has not been already collapsed
    if ismember(keys{i}, [keySets{:}]), continue, end
    
    % find column index
    idx = find(strcmp(dataVars, keys{i}));
    nLevels = length(unique(t_data.(idx)));
    assert(nLevels<height(t_data), 'key %s span the whole table', keys{i})
    
    % find potential matched keys
    for j = setdiff(1:length(dataVars),idx)
        if height(unique(t_data(:,[idx j])))==nLevels 
            %%% unclear how to deal with combined keys and 'annotations'.
            %%% As now it prioritizes the keys given as an input and in the
            %%% order of the input.
            %%%         MH 15/12/18
            %%%% 
            idx = [idx j];
        end
    end
    
    % assign the keys
    keySets{end+1} = dataVars(idx);
    Dimensions{end+1} = unique(t_data(:,idx));
    
end

