function firstdiff = matdiff(a, b)
    ha = size(a, 1);
    hb = size(b, 1);
    h = min(ha, hb);

    l = 1;
    r = h;
    while true
        m = floor((l + r)/2);
        % [l m r]
        if isequal(a(l:m, :), b(l:m, :))
            l = m + 1;
            if l > h; break; end
        else
            if l == m; break; end
            r = m;
        end
    end

    if l > ha; adiff = []; else adiff = [{'a', l} a(l, :)]; end
    if l > hb; bdiff = []; else bdiff = [{'b', l} b(l, :)]; end
    firstdiff = [adiff; bdiff];
end
