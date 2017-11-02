function [reordered_tbl, keyvars, other] = normalize_var_ordering(tbl, varargin)
    narginchk(1, 2);
    if nargin == 1
        args = {tbl};
    else
        args = {tbl 'KeyVars' varargin{1}};
    end

    [~, keyvars] = process_args__({'KeyVars'}, args);
    other = setdiff(varnames(tbl), keyvars, 'stable');
    reordered_tbl = tbl(:, [keyvars other]);
end
