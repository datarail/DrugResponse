function out = togglecr(cellarray)
    dim = getdim__(cellarray);
    assert(~isempty(dim), 'input does not have correct type or size');
    out = num2cell(cat(dim, cellarray{:}), dim);
end

function dim = getdim__(arg)
    dim = [];

    if ~(iscell(arg) && all(cellfun(@iscell, arg))); return; end

    if iscolumn(arg); d = 1;
    elseif isrow(arg); d = 2;
    else return;
    end

    if numel(unique(cellfun(@(rc) size(rc, d), arg))) ~= 1; return; end

    dim = d;
end
