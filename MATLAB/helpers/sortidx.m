function [order, sortedvalues] = sortidx(values, varargin)
% [order, sortedvalues] = sortidx(values, varargin)
%   revert order of output than sort(); same inputs
%

[sortedvalues, order] = sort(values, varargin{:});
