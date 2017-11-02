function [vi vn] = parse_keyvars_arg(tbl, keyvars)

    vns = tbl.Properties.VariableNames;
    function i = idx(v)
        if isnumeric(v)
            i = v;
        else
            [~, i] = ismember(v, vns);
        end
    end

    if ~iscell(keyvars)
        keyvars = {keyvars};
    end

    vi = cellfun(@idx, keyvars);
    unk = keyvars(vi == 0);
    if numel(unk) > 0
        error(['Unrecognized variable(s): ' ...
               strjoin(unique(unk, 'stable'))]);
    end
    if nargout > 1
        vn = vns(vi);
    end
end

