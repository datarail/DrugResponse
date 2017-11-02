function out = tomaxdims(cellarray)
%TOMAXDIMS convert nested cell arrays into a simple cell array.
    if ~iscell(cellarray)
        error(['Input #1 expected to be a cell array, ' ...
               'was %s instead'], class(cellarray));
    end
    out = tomaxdims_(cellarray, 1);
end

function out = tomaxdims_(cellarray, d)
    tf = reshape(cellfun(@iscell, cellarray), 1, []);
    if ~any(tf)
        out = cellarray;
        return;
    end
    assert(all(tf));
    if numel(cellarray) == 1
        out = tomaxdims_(cellarray{1}, d);
        return;
    end
    tmp = arraymap(@(i) tomaxdims_(sqz1_(hslice(cellarray, 1, i)), d+1), ...
                   (1:size(cellarray, 1)).');

    out = cat(d, tmp{:});
end

function out = sqz1_(a)
    sz = size(a);
    assert(sz(1) == 1);
    out = reshape(a, [sz(2:end) 1]);
end

%%
function out = size_(cellarray, d)
    tf = reshape(cellfun(@iscell, cellarray), 1, []);
    if ~any(tf)
        out = size(cellarray);
        return;
    end
    assert(all(tf));
    if numel(cellarray) == 1
        out = size_(cellarray{1}, d);
        return;
    end
    tmp = arraymap(@(i) size_(sqz1_(hslice(cellarray, 1, i)), d+1), ...
                   (1:size(cellarray, 1)).');

    w = max(cellfun(@numel, tmp));
    szs = cell2mat(cellmap(@(s) [s ones(1, w - numel(s))], tmp));
    out = szs(1, :);
    out(d) = sum(szs(:, d));
end
