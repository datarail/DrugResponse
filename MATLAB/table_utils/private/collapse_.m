function out = collapse_(tbl, ag, gv, iv, irreg)
    nargoutchk(0, 3);
    m = numel(iv);
    if m == 0

% $$$     t        = table_([2 3;
% $$$                        1 3;
% $$$                        1 1;
% $$$                        2 1;
% $$$                        1 1;
% $$$                        1 2;
% $$$                        2 2});
% $$$     expected = table_([2 3;
% $$$                        1 3;
% $$$                        1 1;
% $$$                        2 1;
% $$$                        1 2;
% $$$                        2 2});
% $$$     actual = collapse(t, {}, 'KeyVars', {1, 2});

        out = unique(tbl(:, gv), 'rows', 'stable');
    else
        fn = cellmap(@wrapfn, ag);
        try
            out = varfun_(fn, tbl, gv, iv, irreg);
        catch e
            if strfind(e.message, '::DR20:')
                parts = strsplit(e.message, '::');
                exc_id = parts{2};
                exc_msg = parts{3};
                error(exc_id, exc_msg);
            end
            rethrow(e);
        end

% $$$     t        = table_({6 0;
% $$$                        6 0;
% $$$                        4 0;
% $$$                        5 0;
% $$$                        6 0;
% $$$                        5 0});
% $$$     expected = table_({6 3;
% $$$                        4 1;
% $$$                        5 2});
% $$$     actual = collapse(t, @numel, 'KeyVars', {1});
    end
end

%-------------------------------------------------------------------------------
function wrapped = wrapfn(fn)
    function out = wrapped_(x)
        out = fn(x);
        if ~isrow(out)

% $$$     t        = table_({0 1;
% $$$                        0 2});
% $$$     verifyError(testCase, @() collapse(t, @(x) x, 'KeyVars', {1}), ...
% $$$                 'DR20:collapse:ReturnedValueIsNotRowVector');

            error(['::DR20:collapse:ReturnedValueIsNotRowVector::' ...
                   'Output is not a row vector (height = %d)'], ...
                  size(out, 1));
        end
    end
    wrapped = @wrapped_;
end

%-------------------------------------------------------------------------------
function out = unsorter_(t)
    [~, idx] = sortrows(t);
    out(idx, :) = tocol(1:numel(idx));
end

%-------------------------------------------------------------------------------
function out = varfun_(fns, tbl, gvns, vvns, irreg)
    gvidx = dr.vidxs(tbl, gvns);

    nvvs = numel(vvns);
    out_data = cell(1, nvvs);

    [grouped_tbl, starts, ends] = group_rows_(tbl, gvns);

    ngroups = numel(starts);
    nrows = height(grouped_tbl);
    assert(nrows == height(tbl));

%     assert(isequal(tbl(origglocs, :), grouped_tbl(glocs, :)));
%     orig_in_data = tabledata(tbl);
%     for j = 1:size(in_data, 2)
%         assert(isequal(in_data{1, j}(glocs, :), ...
%                        orig_in_data{1, j}(origglocs, :)));
%     end

    out_col = cell(ngroups, 1);
    in_data = tabledata(grouped_tbl);

    vvidx = dr.vidxs(grouped_tbl, vvns);
    for j = 1:nvvs
        vvn = vvns{j};
        vvj = vvidx(j);

        in_col = in_data{vvj};
        fn = fns{j};
        fn_name = func2str(fn);
        on_err = @fn_failed;

        nd = ndims(in_col);
        if nd <= 2
            % common special case
            getarg_ = @(i) in_col(starts(i):ends(i), :);
        else
            template = repmat({':'}, nd);
            getarg_ = @(i) hslice_(in_col, 1, starts(i):ends(i), template);
        end

        for i = 1:ngroups

% $$$     t        = table_([6;
% $$$                        6;
% $$$                        4;
% $$$                        5;
% $$$                        6;
% $$$                        5});
% $$$     to_12345 = @(a) catc(1, arraymap(@(i) i * ones(1:5), a));
% $$$     t.(2) = to_12345(t.(1));
% $$$
% $$$     expected = table_([6 360;
% $$$                        4 120;
% $$$                        5 240});
% $$$
% $$$     actual = collapse(t, @numel, 'KeyVars', {1});

            arg = getarg_(i);
            try out_col{i} = fn(arg);
            catch e
                s = mk_err_struct_(e, fn_name, vvj, vvn, i);
                out_col{i} = on_err(s, arg);
            end
        end

        if irreg{j}
            out_data{j} = out_col;
        else
            try
                out_data{j} = vertcat(out_col{:});
            catch e
                error(message('MATLAB:table:varfun:VertcatFailed', fn_name, ...
                              vvn, e.message));
            end
        end
    end

    assert(nvvs == 0 || all(cellfun(@(x) 1 == size(x, 1), out_col)));

    try
        out = [grouped_tbl(starts, gvidx) ...
               table(out_data{:}, 'VariableNames', vvns)];
    catch e
        e_id = e.identifier;
        if strcmp(e_id, 'MATLAB:table:parseArgs:WrongNumberArgs') || ...
           strcmp(e_id, 'MATLAB:table:parseArgs:BadParamName')
            error(message('MATLAB:table:varfun:CharRowFunOutput'));
        end
        rethrow(e);
    end

end % function b = varfun_(fun, a, varargin)

%-------------------------------------------------------------------------------
function [varargout] = fn_failed(s, varargin)
    throw(MException('MATLAB:table:varfun:FunFailedGrouped', ...
                 ['Applying the function ''%s'' to group %d in ' ...
                 'the variable ''%s'' generated the following ' ...
                 'error:\n\n%s\n'], s.fname, s.group, s.name, s.message));
end

%-------------------------------------------------------------------------------
function s = mk_err_struct_(exc, fname, vvidx, vvname, gidx)
    s = struct('identifier', exc.identifier, 'message', exc.message, ...
               'fname', fname, 'index', vvidx, 'name', vvname, 'group', gidx);
end

%-------------------------------------------------------------------------------
function [glocs, group] = gidx_(ktbl)
    [~, glocs, group] = unique(hashable_(ktbl), 'sorted');
end

%-------------------------------------------------------------------------------
function [glocs, group] = table2gidx_(nrows, kcols, kvns)
    w = numel(kvns);
    if w == 0
        group = ones(nrows, 1);
        glocs = ones(min(nrows, 1), 1);
    elseif w == 1
        [group, glocs] = grp2idx(kcols{1}, kvns{1});
    else
        groups = zeros(nrows, w);
%        for j = 1:w, groups(:, j) = grp2idx(kcols{j}, kvns{j}); end
        for j = 1:w
            groups(:, j) = grp2idx(kcols{j}, kvns{j});
            1;
        end
        assert(~any(isnan(groups(:))));
        [~, glocs, group] = unique(groups, 'rows', 'sorted');
    end
end

%-------------------------------------------------------------------------------
function [gidx, gloc] = grp2idx(var, vname)
    if ischar(var)
        if isempty(var), var = repmat({''}, size(var, 1), 1);
        else var = cellstr(var); end
    end

    %if ~iscolumn(var), error(message('MATLAB:table:GroupingVarNotColumn')); end
    assert(iscolumn(var));
    assert(iscategorical(var) || isnumeric(var) || islogical(var) || ...
           iscell(var));

    if iscategorical(var)
        [glevels, gloc, gidx] = unique(var, 'sorted');
        assert(isempty(glevels) || ~isundefined(glevels(end)));
    else
        try
            [glevels, gloc, gidx] = unique(var, 'sorted');
        catch e
            m = message('MATLAB:table:VarUniqueMethodFailed', vname);
            throwAsCaller(addCause(MException(m.Identifier, getString(m)), e));
        end
        if numel(gidx) ~= numel(var)
            m = message('MATLAB:table:VarUniqueMethodFailedNumRows', vname);
            throwAsCaller(MException(m.Identifier, getString(m)));
        end

        if isempty(glevels), return; end
        if isnumeric(var) || islogical(var)
            assert(~isnan(glevels(end)));
        else
            assert(~isempty(glevels{1}));
        end
    end
    assert(~any(isnan(gidx(:))));
end

