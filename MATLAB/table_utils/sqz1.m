function out = sqz1(a)
    sz = size(a);
    if sz(1) == 1 && numel(sz) > 2
        out = reshape(a, sz(2:end));
    else
        out = a;
    end
end
