function out = perm(x, varargin)
% Y = PERM(X, [CYCLE_1 [, CYCLE_2, ...]])
%
% Y is the result of sequentially permuting the elements of X by the
% cyclic permutations specified by CYCLE_1, CYCLE_2, etc.
%
% >> perm(1:6, [2 6 4], [1 5 3])
% ans =
%      3     4     5     6     1     2

    if nargin == 1
        out = x;
        return;
    end

    out = perm(perm1_(x, varargin{1}), varargin{2:end});
end

function out = perm1_(x, c)
% >> perm1_(1:6, [2 6 4])
% ans =
%      1     4     3     6     5     2
    d = circshift(c, 1, maxdim(c));
    out = x;
    out(c) = out(d);
end
