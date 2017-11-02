function out = transposend( ndarray, varargin )
    narginchk(1, 2);
    p = fliplr(1:ndims(ndarray));
    if nargin > 1 && varargin{1}
        p = circshift(p, [0 -1]);
    end
    out = permute(ndarray, p);
end
