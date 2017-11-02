function tbl = cols2str(tbl)
    for j = 1:width(tbl)
        if iscategorical(tbl.(j))
            tbl.(j) = cellstr(tbl.(j));
        elseif iscell(tbl.(j))
            tbl.(j) = cellmap(@maybe_stringify_, tbl.(j));
        end
    end
end

function c = maybe_stringify_(c)
    if iscell(c)
        c = CStr2String(cellstr(c), '|', 'noTrail');
    end
end
