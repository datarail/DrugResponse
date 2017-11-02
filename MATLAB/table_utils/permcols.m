function out = permcols(tbl, varargin)
% U = PERMCOLS(T, [CYCLE_1 [, CYCLE_2, ...]])
%
% U is the result of sequentially permuting the columns of T by the
% cyclic permutations specified by CYCLE_1, CYCLE_2, etc.

    out = tbl(:, perm(1:width(tbl), varargin{:}));
end
