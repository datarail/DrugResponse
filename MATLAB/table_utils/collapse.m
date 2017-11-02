function out = collapse(tbl, aggrs, varargin)
%COLLAPSE consolidate groups of rows into single rows.
%
%     T2 = COLLAPSE(T1, AGGRS)
%
%     T2 = COLLAPSE(T1, AGGRS, 'PARAM1',val1, 'PARAM2',val2, ...) allows
%     you to specify optional parameter name/value pairs to control
%     COLLAPSE's behavior.  These parameters are listed below.  The
%     documentation of the first two is identical to those of the
%     same-named parameters for the function TABLE_TO_NDARRAY.
%
%         'KeyVars'
%         'ValVars'
%       'Irregular' - Should be either a cell array of boolean values
%                     having as many elements as there are value variables,
%                     or a single boolean variable.  (The second form is a
%                     shorthand for a value in the first form where all the
%                     boolean entries have the value specified.)  The
%                     boolean value corresponding to each value variable
%                     determines how the collapsed values for that variable
%                     will be stacked vertically to form a column of the
%                     returned table.  When this value is TRUE, the entries
%                     will be treated as irregular, and stacked as a
%                     vertical cell array having variably-sized entries.
%                     Otherwise, the entries will be expected to have
%                     uniform dimensions orthogonal to their first
%                     dimension, and therefore suitable for stacking with
%                     VERTCAT.
%
%                     DEFAULT: FALSE.

    narginchk(2, 6);
    args = [{tbl} varargin];
    if ~isempty(aggrs), args = [args {'Aggrs' aggrs}]; end
    [tbl, kns, vns, aggrs, irreg, ~] = ...
        process_args__({'KeyVars' 'ValVars' 'Aggrs' 'Irregular'}, args);
    try
        out = collapse_(tbl, aggrs, kns, vns, irreg);
    catch e
        ds = dbstack('-completenames');
        if isequal(ds(end).name, 'runtests') && ...
           ~any(cellfun(@(s) ...
               ~isempty(strfind(s, 'throwsExpectedException')), {ds.name}.'))
            dbstack;
            rethrow(e);
        end
        exc_id = regexprep(e.identifier, 'MATLAB:.*:', 'DR20:collapse:');
        throw(MException(exc_id, e.message));
    end
end
