function out = table_to_cellstr(tbl)
    out = roc2cor(tabledata(cols2str(tbl)));
end
