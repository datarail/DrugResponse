% simple wrapper for SELECT, to be called when
% 1. SEQ is a cell array, AND
% 2. FN expects its arguments to be passed as a "comma-separated
%    list", rather than as a single cell array
function out = selectc(fn, seq, varargin)
    assert(iscell(seq));
    out = select(x_(fn), seq, varargin{:});
end

function g = x_(f)
    g = @(c) f(c{:});
end