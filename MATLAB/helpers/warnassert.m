function warnassert(condition, varargin)
% warnassert(condition, varargin)
%   through a warning line (using function warnprintf) if condition fails;
%   it works as assert but do not break
%
if ~condition
    warnprintf(varargin{:})
end
