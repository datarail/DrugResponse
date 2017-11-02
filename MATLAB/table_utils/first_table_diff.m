function out = first_table_diff(a, b)
    out = matdiff(table_to_cellstr(a), ...
                  table_to_cellstr(b));
end
