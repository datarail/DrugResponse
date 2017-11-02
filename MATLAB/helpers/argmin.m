function [idx, value] = argmin(varargin)
% [idx, value] = argmin(varargin)
%   invert the output of the min function

[value, idx] = min(varargin{:});
