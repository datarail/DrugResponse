function tbl = ndarray_to_table(ndarray, labels, varargin)
%NDARRAY_TO_TABLE convert table to nd-array.
%     T = NDARRAY_TO_TABLE(A, L) converts nd-array A to table T,
%     having columns as specified by L.
%
%     T = NDARRAY_TO_TABLE(A, L, OUTER) like the previous form, but if
%     OUTER is true, assume that the dimensions of A are ordered as
%     described in the documentation for the OUTER parameter of the
%     function TABLE_TO_NDARRAY.

    narginchk(2, 3);
    outer = nargin > 2 && varargin{1};

    labels = cellmap(@(l) unique(l, 'rows', 'stable'), labels);

    if outer, labels = circshift(labels, [0 1]); end
    vvs = labels{1};
    labels = labels(2:end);

    n = ndims(ndarray);
    check_dims_(n, outer, getnvvs_(vvs), numel(labels));

    if ischar(vvs)
        vvs = {vvs};
    else
        if istable(vvs), vvs = cellstr(vvs.(1)); end
        vvs = reshape(vvs, 1, []);
    end

    tidx = cartesian_product(cellmap(@(t) tocol(1:height(t)), labels));
    keys = catc(2, arraymap(@(i) labels{i}(tidx{i}, :), 1:numel(labels)));

    m = numel(ndarray);

    p = n:-1:1;
    if numel(vvs) > 1
        m = m/numel(vvs);
        if outer, p = circshift(p, [0 -1]); end
    end

    values = num2cell(reshape(permute(ndarray, p), m, []), 1);
    data = cat(2, keys, table(values{:}));
    tbl = make_table(data, keys.Properties.VariableNames, vvs);
end

function check_dims_(ndims_, outer, nvalvars, nkeys)
    % if OUTER is true and the number of valvars is 1, then the number
    % of dimensions should equal the number of keys; otherwise, the
    % number of dimensions should be 1 greater than the number of keys.
    if (outer && nvalvars == 1), d = 0; else d = 1; end
    if ~(ndims_ == d + nkeys)
        error('DR20:ndarray_to_table:InconsistentDimensions', ...
              'Dimensions of ndarray not consistent with label dimensions');
    end
end

function nvvs = getnvvs_(vvs)
    if ischar(vvs), nvvs = 1;
    elseif istable(vvs), nvvs = height(vvs);
    else nvvs = length_(vvs);
    end
end
