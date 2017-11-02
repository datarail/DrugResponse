function out = GET( x, varargin )
%GET functional replacement for parenthesized index expressions.

    out = x(varargin{:});
end

