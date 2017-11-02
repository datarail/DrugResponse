function [n, d] = length_( x, varargin )
%LENGTH_ return length along first non-trivial dimension.
%
%   If the length along some dimension is 1, this dimension is
%   said to be "trivial."
%
%   N = LENGTH_(X) assigns to N the length of the first non-trivial
%   dimension of X.  LENGTH_(X) returns 0 if X is empty, and 1 if all the
%   dimensions of X are trivial.
%
%   If X is a table, assigns HEIGHT(X) to N.
%
%   [N, D] = LENGTH_(X) in addition assigns to D the number of the first
%   non-trivial dimension of X.  This number will be 1 if X is empty or all
%   its dimensions are trivial.
%
%   If X is a table, assigns 1 to D.
%
%   LENGTH_(X, D) is almost equivalent to SIZE(X, D), but it will fail
%   whenever X is a table and D is not 1.  This option serves as a way to
%   specify the length-dimension that is not susceptible to the ambiguities
%   inherent in determining the length dimension from the X object.

    narginchk(1, 2);
    if nargin > 1, d_ = varargin{1}; end

    if istable(x)
        n = height(x);
        if nargin > 1 && d_ ~= 1
          error('DR20:length_:InvalidLengthDimensionForTable', ...
                  ['Specified length dimension (%d) is not valid for ' ...
                   'a table'], d_);
        end
        d = 1;
    else
        sz = size(x);
        if nargin > 1, d = d_;
        else
            if numel(x) > 0
                d = find(sz > 1, 1);
            else
                d = find(sz < 1, 1);
            end
            if isempty(d), d = 1; end
        end
        n = sz(d);
    end
end
