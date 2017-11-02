function out = grouped_select(tbl, selector, varargin)
%GROUPED_SELECT consolidate groups of rows into single rows.
%
%     T2 = GROUPED_SELECT(T1, SELECTOR)
%
%     T2 = GROUPED_SELECT(T1, SELECTOR, 'PARAM1',val1, 'PARAM2',val2, ...) allows
%     you to specify optional parameter name/value pairs to control
%     GROUPED_SELECT's behavior.  These parameters are listed below.  ...
%     documentation of the first two is identical to those of the
%     same-named parameters for the function TABLE_TO_NDARRAY.
%
%         'KeyVars'
%         'ValVar'

    narginchk(2, 6);
    args = [{tbl} varargin];
    [tbl, kns, vns, ~, ~, ~, parser] = ...
        process_args__({'KeyVars' 'ValVars'}, args);

    if ismember('ValVars', p.UsingDefaults); then vns = {}; end

    out = grouped_select_(tbl, selector, kns, vns);
end
