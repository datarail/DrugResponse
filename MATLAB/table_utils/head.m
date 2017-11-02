function out = head(seq, varargin)
% out = head(seq, varargin)
    narginchk(1, 2);
    if nargin > 1
        n = varargin{1};
        if n < 0
            n = numel(seq) + n;
        end
        assert(n > 0);
    else
        n = 10;
    end
    out = take(seq, n);
end
