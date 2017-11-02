function [idx, value] = argmax(varargin)
% [idx, value] = argmax(varargin)
%   invert the output of the max function

[value, idx] = max(varargin{:});
