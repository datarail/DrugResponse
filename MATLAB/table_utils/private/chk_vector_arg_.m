function [] = chk_vector_arg_(x, nargin_, varargin_)
    assert(nargin == 3);
    if ~(isvector_(x) || (nargin_ > 1 && varargin_{1}))
        error(['DR20:' caller ':ArgNotVector'], 'Argument is not a vector');
    end
end
