function out = cellmap(fun, c)
%CELLMAP(FUN, C) is a shorthand for CELLFUN(FUN, C, 'UniformOutput', false)
    out = cellfun(fun, c, 'UniformOutput', false);
end
