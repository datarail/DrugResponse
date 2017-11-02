function out = group_counts(tbl, varargin)
%GROUP_COUNTS count numbers of rows in groupings.
%
%     T2 = GROUP_COUNTS(T1 [, KEYVARS])
%
%     Returns a table whose columns correspond to KEYVARS plus an
%     additional column holding the numbers of rows in T1 having each
%     available combination of values of the variables in KEYVARS.
%
%     If KEYVARS is omitted, the default is computed as described in
%     the documentation for the KEYVARS parameter of the function
%     TABLE_TO_NDARRAY.
%

    narginchk(1, 2);
    if nargin > 1, params = {'KeyVars', varargin{1}}; else params = {}; end

    [tbl, kvs, ~, ~, ~, ~] = process_args__({'KeyVars'}, [{tbl} params]);

    kis = dr.vidxs(tbl, kvs);
    t = tbl(:, kis);
    vv = genvarname_(kvs, 'counts');
    t.(vv) = ones(height(t), 1);

    out = collapse_(t, {@numel}, kvs, {vv}, {false});
end
