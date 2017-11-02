function tbl = dropcols(tbl, dropset)
%DROPCOLS Remove some columns.
%   B = DROPCOLS(A, COLSTODROP) produces a new table B
%   with the same columns as those in A (in the same order), except
%   for those that are listed in COLSTODROP.
%
%   The class (table) of the first argument determines
%   the class of the result.
%
%   The second argument may be a vector of integers, representing
%   column indices, or a cell array of column names and/or column
%   indices.  The ordering of these column designators is not
%   significant.  Repeated designators are ignored.  Unrecognized
%   designators trigger an exception.

    w = size(tbl, 2);
    tbl = tbl(:, setdiff(1:w, dr.vidxs(tbl, dropset)));
end
