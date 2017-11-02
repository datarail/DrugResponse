function d = tabledata(t)
%TABLEDATA Return contents of table as a cell array of columns.
%   B = TABLEDATA(A) is shorthand for
%
%   B = arrayfun(@(i) A.(i), 1:width(A), 'UniformOutput', false);

    d = arraymap(@(i) t.(i), 1:width(t));
end
