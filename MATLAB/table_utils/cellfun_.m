function output = cellfun_(fhandle, varargin)
% output = cellfun_(fhandle, varargin)
%   same use as the cellfun function with non-uniform output

output = cellfun(fhandle, varargin{:}, 'uniformoutput', false);
