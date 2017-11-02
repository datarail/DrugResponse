function out = catc(dim, cellarray)
%CATC Concatenate cell arrays.
%   B = CATC(DIM, CELLARRAY) is equivalent to B = CAT(DIM, CELLARRAY{:}).

    if iscell(cellarray)
        out = cat(dim, cellarray{:});
    elseif isscalar_(cellarray)
        out = cellarray;
    else
        error('DR20:catc:UnsupportedClass', 'Unsupported class: %s', ...
               class(cellarray));
    end
end
