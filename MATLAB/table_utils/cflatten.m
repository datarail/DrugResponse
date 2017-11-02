function out = cflatten(c)
    if iscell(c)
        out = catc(2, cellmap(@cflatten, reshape(c, 1, [])));
    else
        out = {c};
    end
end
