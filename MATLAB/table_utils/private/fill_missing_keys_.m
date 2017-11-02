function P = fill_missing_keys_(tbl, kns)

    levels = cellmap(@(n) unique(tbl.(n), 'stable'), kns);

    max_unique = prod(cellfun(@numel, levels));
    if height(unique(tbl(:, kns))) == max_unique
        P = tbl;
        return;
    end
    clear('max_unique');

    s = cartesian_product_table(levels, kns);

    %%% unstable ordering
    % P = outerjoin(s, tbl, 'MergeKeys', true);

    %%% stable ordering
    [P, itbl, chk] = outerjoin(tbl, s, 'MergeKeys', true);
    assert(all(chk > 0));
    clear('s', 'chk');

    % the non-zero entries in itbl correspond to tbl's original ordering;
    % the code below replaces the zero entries in itbl with consequtive
    % indices following the greatest non-zero index in itbl; this will
    % result a final ordering in which the original rows, in their
    % original order, are followed by the new rows.

    missing = find(itbl == 0); assert(~isempty(missing));
    h = height(P);
    g = h - numel(missing) + 1;
    assert(all(itbl < g));
    itbl(missing) = (g:h).';
    idx(itbl, 1) = (1:h).';
    P = P(idx, :);
end
