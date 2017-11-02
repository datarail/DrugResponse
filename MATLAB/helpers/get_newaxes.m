function ahandle = get_newaxes(position, holded, varargin)
% ahandle = get_newaxes(position, holded, varargin)
%   create a new axis with the given position

if nargout==0
    axes('position', position, varargin{:});
else
    ahandle = axes('position', position, varargin{:});
end
if exist('holded','var') && ~isempty(holded) && holded
    hold on
end
