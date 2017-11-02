function [grouped_tbl, starts, ends, sort_perm, sort_iperm] = group_rows(tbl, varargin)
%GROUPROWS reorder TBL's rows to group them by key-variable values.
%
%     T2 = GROUPROWS(T1, KEYVARS)
%     T2 = GROUPROWS(T1)
%
% The result is similar to that of SORTROWS(T1, KEYCOLS), but
% GROUPROWS orders the groups in the order of first appearance of
% distinct value combinations in the key variables.

    narginchk(1, 2);
    if nargin == 1
        args = {tbl};
    else
        args = {tbl 'KeyVars' varargin{1}};
    end

    [tbl, kns] = process_args__({'KeyVars'}, args);
    [grouped_tbl, starts, ends, sort_perm, sort_iperm] = group_rows_(tbl, kns);
end
