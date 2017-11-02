function arr = make_test_ndarray( sz, varargin )
    narginchk(1, 2);
    sz = reshape(sz, 1, []);
    if ~all(sz > 0)
        error('first argument contains zero dim(s)');
    end
    n = numel(sz);
    if nargin == 2 && varargin{1}
        % arr = dr.mkbox(sz, true, true);
        % next:
        arr = dr.mkbox(sz, false, false);
    else
        % arr = dr.mkslab(sz, true, true);
        % next:
        arr = dr.mkslab(sz, false, false);
    end
end
