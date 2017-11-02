function out = table2mat_( tbl, vns )
    cols = cellmap(@(i) tbl.(i), vns);
    out = cat(2, cols{:});
end

