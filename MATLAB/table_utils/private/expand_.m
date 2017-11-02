function out = expand_(n)
    out = cast(arrayfun(@str2double, int2str(n)), class(n));
end
