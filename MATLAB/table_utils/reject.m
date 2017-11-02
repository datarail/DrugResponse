function out = reject( fn, seq )
    [out, sz, dn, idx] = filter_0_(fn, seq);
    out = filter_1_(out, sz, dn, ~idx);
end

