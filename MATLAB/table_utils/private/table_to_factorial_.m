function out = table_to_factorial_(tbl, kns, vns, aggrs)
    irreg = repmat({false}, 1, numel(vns));
    out = fill_missing_keys_(collapse_(tbl, aggrs, kns, vns, irreg), kns);
end
