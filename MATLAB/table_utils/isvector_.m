function tf = isvector_(x)
    sz = size(x);
    %tf = numel(sz(sz > 1)) <= 1 && ~any(sz == 0);
    tf = numel(sz(sz > 1)) <= 1 && numel(sz(sz == 1)) >= 1;
end
