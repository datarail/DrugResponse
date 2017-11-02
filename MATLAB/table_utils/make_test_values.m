function out = make_test_values( sz, varargin )
%MAKE_TEST_VALUES produce array of test values.

    narginchk(1, 3);
    sz = reshape(sz, 1, []);
    if ~all(sz > 0)
        error('first argument contains zero dim(s)');
    end
    expanded = nargin > 1 && varargin{1};
    colmajor = nargin > 2 && varargin{2};

    if colmajor
        out = fliplr(mtv_(fliplr(sz), numel(sz)));
    else
        out = mtv_(sz, numel(sz));
    end
    if ~expanded
        out = arrayfun(@(i) contract_(out(i, :)), (1:size(out, 1)).');
    end
end

function out = mtv_( sz, n )
    if n == 1
        out = (1:sz(1)).';
    else
        base = mtv_(sz(2:end), n - 1);
        n = size(base, 1);
        tmp = arraymap(@(i) [repmat(i, n, 1) base], 1:sz(1));
        out = vertcat(tmp{:});
    end
end
