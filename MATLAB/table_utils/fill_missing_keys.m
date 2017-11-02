function P = fill_missing_keys(varargin)
%FILL_MISSING_KEYS(TBL, KEYVARS) produces a new table from the table TBL by
%    inserting new rows for all the combinations of possible values of
%    the KEYVARS variables that are not present in TBL.  The values
%    for the remaining variables of TBL are those produced by OUTERJOIN.
%
%FILL_MISSING_KEYS(TBL) is equivalent to a call to
%    FILL_MISSING_KEYS(TBL, KEYVARS) with KEYVARS set to hold the
%    variables K of TBL for which ISCATEGORICAL(TBL.K) is true.

    % [tbl, kis, ~, ~, ~] = process_args__({'KeyVars'}, varargin);
    % kns = dr.vns(tbl, kis);
    [tbl, kns, ~, ~, ~] = process_args__({'KeyVars'}, varargin);
    P = fill_missing_keys_(tbl, kns);
end
