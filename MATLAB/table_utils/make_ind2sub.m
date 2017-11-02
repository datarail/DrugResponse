function fn = make_ind2sub( sz, varargin )
    narginchk(1, 2);
    n = numel(sz);

    function out = i2s_0(i)
        [idx{1:n}] = ind2sub(sz, i);
        out = cell2mat(idx);
    end
    if nargin > 1
        cls = varargin{:};
    else
        cls = '';
    end
    function out = i2s_1(i)
        [idx{1:n}] = ind2sub(sz, i);
        out = cast(cell2mat(idx), cls);
    end

    if nargin == 1
        fn = @i2s_0;
    else
        fn = @i2s_1;
    end

end

