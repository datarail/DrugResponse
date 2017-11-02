function out = slice1_( cellarray, widths, varargin )
%SLICE1_ carve array into subarrays.
%    ARRAY is either an n-dimensional array, or a regular cell array
%    containing only non-cell elements of a uniform type (as would be
%    produced, for example, by applying NUM2CELL to an n-dimensional
%    array).
%
%    WIDTHS is a array such that SIZE(WIDTHS) is [1
%    NDIMS(ARRAY)], and whose I-th entry is a row vector of
%    positive integers whose sum does not exceed SIZE(ARRAY, I).
%
%    B = slice1_(A, WIDTHS) sets B to a array such that SIZE(B)
%    is equal to CELLFUN(@NUMEL, WIDTHS), and that SIZE(B{I_1, ... ,
%    I_n}) equals [WIDTHS{1}(I_1) ... WIDTHS{n}(I_n)].
%
%    C = slice1_(A, WIDTHS, true) sets C to a array such that SIZE(C)
%    is equal to NUMEL(WIDTHS{1}), and that SIZE(C{I_1}{...}{I_n})
%    equals [WIDTHS{1}(I_1) ... WIDTHS{n}(I_n)].
%
%    The essential difference between the two forms shown above is that,
%    for the first one, the slicing is done simultaneously across all the
%    dimensions of A, whereas for the second one, the slicing is done
%    serially (the result is obtained by splitting A along its first
%    dimension, recursively applying the same serial procedure to the
%    resulting pieces, and concatenating the results.
%
%    Example
%    If A, B, and C are initialized as follows:
%
%    A = [11 12 13 14 15 16; ...
%         21 22 23 24 25 26; ...
%         31 32 33 34 35 36; ...
%         41 42 43 44 45 46; ...
%         51 52 53 54 55 56; ...
%         61 62 63 64 65 66];
%
%    B = slice1_(A, {[1 2 3] [1 2 3]});
%    C = slice1_(A, {[1 2 3] [1 2 3]}, true);
%
%    ...then each B{i, j} and each C{i}{j} is an i-by-j cell array.  E.g.
%
%    >> cell2mat(B{2, 3})
%    ans =
%        24    25    26
%        34    35    36
%
%    >> cell2mat(C{2}{3})
%    ans =
%        24    25    26
%        34    35    36
%
%    In fact, B will be a 3-by-3 cell array:
%
%    B =
%        {1x1 cell}    {1x2 cell}    {1x3 cell}
%        {2x1 cell}    {2x2 cell}    {2x3 cell}
%        {3x1 cell}    {3x2 cell}    {3x3 cell}
%
%    ...and C will be a 3-by-1 cell array:
%
%    C =
%        {3x1 cell}
%        {3x1 cell}
%        {3x1 cell}

    narginchk(2, 3);
    serial = nargin > 2 && varargin{1};

    [cellarray, idxs] = process_args_(cellarray, widths);

    if serial
        out = serial_slice_(cellarray, idxs, 1);
    else
        out = parallel_slice_(cellarray, idxs);
    end
end

function out = serial_slice_(cellarray, idxs, d)
    if isempty(idxs)
        out = cellarray;
    else
        out = cellmap(@(c) serial_slice_(c, idxs(2:end), d+1), ...
                      cellmap(@(ix) hslice(cellarray, d, ix), ...
                      idxs{1}.'));
    end
end

function out = parallel_slice_(cellarray, idxs)
    sz = cellfun(@numel, idxs);
    nd = numel(sz);
    out = cell(sz);
    cp = cartesian_product(idxs, true);
    cp = cat(2, cp{:});
    for i = 1:size(cp, 1)
        subs = cp(i, :);
        [ii{1:nd}] = ind2sub(sz, i);
        out(ii{:}) = {cellarray(subs{:})};
    end
end

function [c, idxs] = process_args_(cellarray, widths)
    if ~iscell(cellarray)
        cellarray = num2cell(cellarray);
    end
    assert(all(reshape(cellfun(@(e) ~iscell(e), cellarray), ...
                       [], 1)), ...
           'argument #1 is not a simple cellarray');
    assert(isvector(widths) && ...
           ndims(cellarray) == numel(widths) && ...
           all(cellfun(@(w) isvector(w) && ~iscell(w), ...
           (widths(:)).')), ...
           'argument #2 is malformed');
    widths = reshape(widths, 1, []);
    assert(all(cellfun(@(w) all(isint(w) & (w > 0)), widths)), ...
           'argument #2 contains invalid widths');
    widths = cellmap(@(w) reshape(w, 1, []), widths);
    ds = size(cellarray) - cellfun(@sum, widths);
    assert(all(ds >= 0), 'some widths exceed size of argument #1');
    if (any(ds > 0))
        widths = cellmap(@(w) w(w ~= 0), ...
                         arraymap(@(i) [widths{i} ds(i)], ...
                                  1:numel(ds)));
    end
    c = cellarray;
    idxs = cellmap(@toidxs_, widths);
end

function idxs = toidxs_(ws)
    assert(~iscell(ws) && isvector(ws) && size(ws, 1) == 1);
    nw = numel(ws);
    idxs = cell(1, nw);
    from = 0;
    for i = 1:nw
        to = from + ws(i);
        assert(to > from);
        idxs{i} = from+1:to;
        from = to;
    end
end
