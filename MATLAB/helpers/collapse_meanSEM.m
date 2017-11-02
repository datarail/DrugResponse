function t_out = collapse_meanSEM(t_in, keys, valvars, pre_trans, applied_fcts)
% t_out = collapse_meanSEM(t_in, keys, valvars, pre_trans, applied_fcts)
%   based on collpase function; by default the @mean and @SEM is perfomed
%   on each value variable (i.e. not key varaible)
%
%   Options:
%       the value variables can be transformed prior to mean/SEM with the
%           function pointed by pre_trans
%       the functions can be changed with the optional input applied_fcts
%

if iscellstr(keys) || ischar(keys)
    keyidx = find(ismember(t_in.Properties.VariableNames, keys));
else
    keyidx = keys;
end
if exist('valvars','var') && ~isempty(valvars) && (iscellstr(valvars) || ischar(valvars))
    validx = find(ismember(t_in.Properties.VariableNames, valvars));
else
    validx = setdiff(1:size(t_in,2), keyidx);
    validx = validx(all(cell2mat(cellfun2(@(x) isnumeric(x) && isscalar(x), ...
        table2cell(t_in(1:min(end,5), validx))))));
end
if exist('pre_trans','var') && ~isempty(pre_trans)
    for i=1:length(validx)
        t_in.(validx(i)) = pre_trans(t_in.(validx(i)));
    end
end

if ~exist('applied_fcts','var')
    applied_fcts = {@mean @SEM};
end

t_out = collapse(t_in(:,[keyidx repmat(validx, 1, length(applied_fcts))]), ...
    reshape(repmat(applied_fcts,length(validx),1),1,[]), ...
    'keyvars', 1:length(keyidx), ...
    'valvars', length(keyidx)+(1:(length(applied_fcts)*length(validx))));

for i=1:length(applied_fcts)
    idx = length(keyidx)+(1:length(validx))+(i-1)*length(validx);
    if i==1
        t_out.Properties.VariableNames(idx) = strcat(t_out.Properties.VariableNames(idx), ...
            ['_' func2str(applied_fcts{i})]);
    else
        t_out.Properties.VariableNames(idx) = strcat(...
            cellfun2(@(x) x(1:(end-2)), t_out.Properties.VariableNames(idx)), ...
            ['_' func2str(applied_fcts{i})]);
    end
end
