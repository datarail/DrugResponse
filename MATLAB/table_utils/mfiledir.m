function out = mfiledir(varargin)
    narginchk(0, 1);
    if nargin == 0
      offset = 0;
    else
      offset = varargin{1};
    end
    out = fileparts(mfilepath(offset + 1));
end
