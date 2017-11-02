function output = cellfun2(fhandle, varargin)
% output = cellfun2(fhandle, varargin)
%   same use as the cellfun function with non-uniform output

output = cellfun(fhandle, varargin{:}, 'uniformoutput', false);
