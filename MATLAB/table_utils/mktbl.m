function out = mktbl(rows, varnames)
    assert_(size(rows, 2) == 1, 'first argument is not a column');
    check_rows_(rows);
    w = get_width_(rows);
    assert_(numel(varnames) == w, 'inconsistent widths');

    varnames = reshape(varnames, 1, []);
    rows = reshape(rows, [], 1);
    cols = arraymap(@(i) ith_col_(rows, i), 1:w);
    out = table(cols{:}, 'VariableNames', varnames);
end

function out = ith_col_(rows, i)
    assert(size(rows, 2) == 1);
    ccol = cellmap(@(r) r{i}, rows);
    if ~isstring_(ccol) && isuniform_(ccol)
        try
            out = catc(1, ccol);
            return;
        catch e
            e
            % ignore
        end
    end
    out = ccol;
end

function out = isuniform_(ccol)
    out = (numel(unique(cellmap(@class, ccol))) < 2);
end

function out = isstring_(ccol)
    out = numel(ccol) > 0 && isstr_(ccol{1});
end

function out = get_width_(rows)
    out = unique(cellfun(@(r) size(r, 2), rows));
    assert_(numel(out) == 1, 'first argument is ragged');
end

function out = check_rows_(rows)
    hs = unique(cellfun(@(r) size(r, 1), rows));
    assert_(numel(hs) == 1 && hs == 1, 'elements are not rows');
end

function assert_(pred, msg)
    assert(pred, 'DR20:failed_assertion', msg)
end
