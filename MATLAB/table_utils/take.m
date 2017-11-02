function out = take(seq, n)
%TAKE Return first N elements from sequence.

    [l, d] = length_(seq);
    if n > l
        out = seq;
    elseif istable(seq)
        out = seq(1:n, :);
    else
        sz = size(seq);
        sz(d) = n;
        out = reshape(seq, l, []);
        out = reshape(out(1:n, :), sz);
    end
end
