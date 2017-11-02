function out = torow(x, varargin)
    narginchk(1, 2);
    chk_vector_arg_(x, nargin, varargin);
    out = reshape(x, 1, []);
end
