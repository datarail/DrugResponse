function [ndarray, labels] = table_to_ndarray(varargin)
%TABLE_TO_NDARRAY convert table to nd-array.
%     A = TABLE_TO_NDARRAY(TBL) converts table TBL to nd-array A.
%
%     [A, L] = TABLE_TO_NDARRAY(TBL) also puts the labels for A into
%     the cell array L.
%
%     A = TABLE_TO_NDARRAY(TBL, 'PARAM1',val1, 'PARAM2',val2, ...) allows
%     you to specify optional parameter name/value pairs to control
%     TABLE_TO_NDARRAY's behavior.  Parameters are:
%
%         'KeyVars' - Variables of TBL to use as key variables.  Only
%                     variables that correspond to 1-dimensional table
%                     columns may be specified as key variables.  The
%                     ordering of the variables in this parameter
%                     affects the ordering of the dimensions of the
%                     resulting NDARRAY.  Key variables may be
%                     specified as a one-dimensional cell array
%                     containing variable names, variable numbers, or
%                     combination thereof, or some one-dimensional
%                     value (e.g. a numeric vector) resolvable to such.
%                     If, after resolving all the specified variables,
%                     duplicates are detected, an exception is raised.
%
%                     DEFAULT: variables of TBL that are not specified
%                     as ValVarS (see below) and for which the ISCATEGORICAL
%                     predicate returns TRUE (ordered according to their
%                     original ordering among TBL's variables).
%
%         'ValVars' - Variables of TBL to use as value variables.  The
%                     ordering of the variables in this parameter
%                     determines the ordering of the corresponding slices
%                     in the resulting NDARRAY.  Value variables may be
%                     specified as a one-dimensional cell array containing
%                     variable names, variable numbers, or combination
%                     thereof, or some one-dimensional value (e.g. a
%                     numeric vector) resolvable to such.  Repeated value
%                     variables are preserved (and will result in repeated
%                     slices in the resulting NDARRAY).
%
%                     DEFAULT: non-key variables of TBL for which the
%                     ISCATEGORICAL predicate returns FALSE, in their
%                     original order of appearance.
%
%           'Aggrs' - Either a function handle or a cell array of function
%                     handles.  This parameter is required if the table
%                     contains multiple rows for some combination of values
%                     of the KEYVARS.  The function(s) specified by this
%                     parameter will be used to aggregate the values of the
%                     corresponding VALVARS for each group of rows having
%                     the same combination of values of the KEYVARS.  If a
%                     function handle is specified, it is used for all the
%                     value variables.  If a cell array of function handles
%                     is specified, it must contain one function handle for
%                     each value variable.  The function(s) specified
%                     should not be sensitive to the ordering of argument
%                     values, and must return values for which ISROW
%                     returns true.
%
%                     DEFAULT: (applicable only if every combination of
%                     values of the KEYVARS is unique) the identity
%                     function.
%
%           'Outer' - A logical value: if TRUE (FALSE), the slices in the
%                     resulting NDARRAY will corresponding to indices of
%                     NDARRAY's last (first) dimension.  Note that if the
%                     number of value variables is 1 and the value of this
%                     parameter is TRUE, the resulting NDARRAY will have K
%                     dimensions (where K is the number of key variables);
%                     otherwise it will have K + 1 dimensions.  The reason
%                     for this is that MATLAB automatically neglects
%                     trailing dimensions of size 1.
%
%                     DEFAULT: FALSE.

    [tbl, kns, vns, aggrs, ~, outer] = ...
        process_args__({'KeyVars' 'ValVars' 'Aggrs' 'Outer'}, varargin);

    ftbl = table_to_factorial_(tbl, kns, vns, aggrs);

    [ii, levels] = tbl_ndarray_ordering_(ftbl, dr.vidxs(ftbl, kns));

    sh = cellfun(@numel, levels);
    assert(prod(sh) == height(ftbl), 'first argument is not factorial');

    nvs = numel(vns);

    vlvls = categorical(vns);
    if outer
        sh = [sh nvs];
        ndarray = reshape(ftbl{ii, vns}, sh);
        levels = [levels {vlvls}];
    else
        sh = [nvs sh];
        ndarray = reshape(ftbl{ii, vns}.', sh);
        levels = [{vlvls} levels];
    end

    if nargout > 1
        if outer
            lns = [kns {'Value'}];
        else
            lns = [{'Value'} kns];
        end
        nls = numel(levels);
        assert(numel(lns) == nls);
        labels = arraymap(@(i) tolabels_(levels{i}, lns{i}), 1:nls);
    end

end


function [ii, levels] = tbl_ndarray_ordering_(tbl, kis)
    m = max(kis);
    assert(0 < min(kis) && m <= width(tbl), ...
        'second argument contains invalid key indices');
    ki2levels = cell(m);
    for ki = kis
        ki2levels{ki} = unique(tbl{:, ki}, 'stable');
    end

    levels = ki2levels(kis);
    % the fliplr below yields a colexicographic ordering wrt
    % the kis;
    svs = arraymap(@(ki) dr.indexof(tbl{:, ki}, ki2levels{ki}), ...
        fliplr(kis));
    [~, ii] = sortrows(cell2mat(svs));
end

function out = tolabels_(lbls, name)
    out = table(lbls(:), 'VariableNames', {name});
end
