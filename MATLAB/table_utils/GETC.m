function out = GETC( x, varargin )
%GETC functional replacement for curly-brace index expressions.

    out = x{varargin{:}};
    %out = builtin('_brace', ca, varargin{:});
    %out = subsref(x, struct('type', '{}', 'subs', {varargin}));
end

