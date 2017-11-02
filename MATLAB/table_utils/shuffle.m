function out = shuffle(x, varargin)
    out = reshape(x(randperm(numel(x))), size(x));
end
