function [varargout] = process_args__(tokeep, orig)

    keep = @(nm) ismember(nm, tokeep);

    assert(~(keep('Aggrs') || keep('Irregular')) || keep('ValVars'));

    p = inputParser;
    p.addRequired('tbl');
    params = {'KeyVars',   []; ...
              'ValVars',   []; ...
              %'Aggrs',     @(x) sum(x, 'native'); ...
              %'Aggrs',     @(x) x(1); ...
              'Aggrs',     @fail_on_repeats; ...
              'Irregular', false; ...
              'Outer',     false};

    for kv = params.'
        if keep(kv{1}), p.addParameter(kv{:}); end
    end

    function [varargout] = unpack_ca_(c)
        [varargout{1:nargout}] = c{:};
    end

    defaults = params(:, 2).';
    [kis, vis, aggrs, irregular, outer] = unpack_ca_(defaults);
    clear('defaults');
    p.parse(orig{:});
    args = p.Results;

    tbl = args.tbl;
    w = width(tbl);
    if keep('KeyVars')
        default_keyvars = ismember('KeyVars', p.UsingDefaults);
        if default_keyvars
            kis = torow(select(@(i) iscategorical(tbl.(i)), 1:w));
        else
            kis = dr.vidxs(tbl, args.KeyVars);
        end
    else
        default_keyvars = true;
    end

    % TODO: implement support for allowing the setting of ValVars to
    % be a function handle that returns a logical scalar.  In this
    % case, only those columns in TBL for which that function returns
    % true (and are not among user-specified keyvars, if any) are
    % treated as valvars.
    % TODO: update documentation accordingly
    %
    %     if isa(args.ValVars, 'function_handle')
    %         try
    %             vis = find(cellfun(args.ValVars, data_(tbl)));
    %         catch e
    %             if strcmp(e.identifier, 'MATLAB:cellfun:NotAScalarOutput')
    %                 error(err_id_('InvalidValVarFunc'), ...
    %                      'Invalid ValVar function')
    %             else
    %                 rethrow(e);
    %             end
    %         end
    %         if default_keyvars
    %             kis = setdiff(kis, vis, 'stable');
    %         else
    %             vis = setdiff(vis, kis, 'stable');
    %         end
    %     end
    %
    % ...
    %
    % function d = data_(t)
    %     d = arraymap(@(i) t.(i), 1:width(t));
    % end

    if keep('ValVars')
        if ismember('ValVars', p.UsingDefaults)
            vis = torow(setdiff(select(@(i) ~iscategorical(tbl.(i)), 1:w), ...
                                kis, 'stable'));
        else
            vis = dr.vidxs(tbl, args.ValVars);

            if default_keyvars
                kis = setdiff(kis, vis, 'stable');
            else
                assert(isempty(intersect(vis, kis, 'stable')), ...
                       'non-empty keyvar-valvar intersection');
            end

            if ~default_keyvars
                assert(isempty(intersect(vis, kis, 'stable')), ...
                       'non-empty keyvar-valvar intersection');
            end
        end
    end

    nvis = numel(vis);

    if keep('Aggrs')
        assert(keep('ValVars'));
        aggrs = args.Aggrs;
        if iscell(aggrs)
            % if make it to here: user specified array of aggrs
            assert(numel(aggrs) == nvis, ...
                   'numbers of valvars and aggregators differ');
            invalid_ag = cellfun(@(f) ~isa(f, 'function_handle'), aggrs);
        else
            invalid_ag = ~isa(aggrs, 'function_handle');
            if ~invalid_ag, aggrs = repmat({aggrs}, 1, nvis); end
        end
        if any(invalid_ag)
            error('InvalidFunction', ...
                  'At least one specified aggregator is not a valid function');
        end
        clear('invalid_ag');
    end

    if keep('Irregular')
       irregular = args.Irregular;
       if ~iscell(irregular), irregular = repmat({irregular}, 1, nvis); end
    end

    if keep('Outer'), outer = args.Outer; end

    kns = dr.vns(tbl, kis);
    vns = dr.vns(tbl, vis);

    missing = select(tooneargfn(@(i) hasmissing_(tbl.(i), i)), kns);
    if ~isempty(missing)
        disp(missing)
        exc = MException(err_id_('MissingValuesInKeyCols'), ...
                         'Some key columns have missing values: %s', ...
                         strjoin(missing, ', '));
        throwAsCaller(exc);
    end

    out = {tbl kns vns aggrs irregular outer p};
    [varargout{1:nargout}] = out{1:nargout};
end

function x = fail_on_repeats(x)
    if size(x, 1) ~= 1
        error(['keyvars do not make up a key: ', ...
               'aggregator(s) must be specified']);
    end
end

function id_ = err_id_(sfx, varargin)
    if nargin > 1, offset = varargin{1}; else offset = 0; end
    id_ = strjoin({'DR20' caller(offset + 1) sfx}, ':');
end

function tf = hasmissing_(var, varname)
    if ischar(var)
        if isempty(var), var = repmat({''}, size(var, 1), 1);
        else var = cellstr(var); end
    end

    if islogical(var)
        assert(~any(isnan(var)));
        tf = false;
    elseif isnumeric(var)
        tf = any(isnan(var));
    elseif iscategorical(var)
        tf = any(isundefined(var));
    elseif iscell(var)
        tf = any(strcmp(var, ''));
    else
        exc = MException(err_id_('UnsupportedKeyVarClass', 1), ...
                         'Class of key variable %s is not supported: %s', ...
                         varname, class(var));
        throwAsCaller(exc);
    end

    tf = isscalar(tf) && tf;
end
