function out = filter_1_( out, sz, dn, idx )
    if any(~idx)
        out = out(idx, :);
        sz(dn) = size(out, 1);
    end
    out = reshape(out, sz);
end
