function out = mfilepath(varargin)
    narginchk(0, 1);
    if nargin == 0
      offset = 0;
    else
      offset = varargin{1};
    end

    tmp = dbstack(1 + offset, '-completenames');
    out = tmp(1).file;
end
