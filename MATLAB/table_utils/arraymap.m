function out = arraymap(fun, b)
%ARRAYMAP(FUN, B) is a shorthand for ARRAYFUN(FUN, B, 'UniformOutput', false)
    out = arrayfun(fun, b, 'UniformOutput', false);
end
