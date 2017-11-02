function out = select( fn, seq, varargin )
    [n, d] = length_(seq, varargin{:});
    idx = arrayfun(@(i) fn(hslice(seq, d, i)), (1:n));
    out = hslice(seq, d, find(idx));
end
