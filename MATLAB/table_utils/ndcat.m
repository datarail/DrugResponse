function ndarray = ndcat(ndas, varargin)
    % NDCAT(NDAS)
    % NDCAT(NDAS, false)
    % -> EXTERNAL indexing

    % NDCAT(NDAS, true)
    % -> INTERNAL indexing
    narginchk(1, 2);
    n = numel(ndas);
    if n == 0
        ndarray = [];
        return;
    end
    s1 = size(ndas{1});
    if n > 1
        chkszs_(s1, ndas(2:end));
    end
    d = numel(s1) + 1;
    ndarray = cat(d, ndas{:});
    if nargin > 1 && varargin{1}
        ndarray = permute(ndarray, circshift(1:d, [0 1]));
    end
end

function chkszs_(sz, ndas)
    assert(all(cellfun(@(a) isequal(sz, size(a)), ndas)));
end
