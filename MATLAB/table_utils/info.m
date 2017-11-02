function [] = info(varargin)
    function s = fmt_(x)
        cls = class(x);
        n = max(2, 12 - numel(cls));
        m = 5;
        fmt = sprintf('%%s %%%dd %%%dd', n, m);
        sz = size(x);
        s = sprintf(fmt, cls, sz(1), sz(2:end));
    end
    records = cellmap(@fmt_, varargin.');
    %disp(char(13));  % trick disp into producing a single newline
    if numel(records) > 1
        nr(records, 'noquote');
    else
        disp(sprintf('\t%s', records{:}));
    end
    %disp(char(13));
end
