function [out, sz, d, idx] = filter_0_( fn, seq )
    [n, d] = length_(seq);

    sz = size(seq);
    out = reshape(seq, n, []);

    if iscell(seq)
        idx = arrayfun(@(i) fn(out{i, :}), (1:n).');
    else
        idx = arrayfun(@(i) fn(out(i, :)), (1:n).');
    end
end
