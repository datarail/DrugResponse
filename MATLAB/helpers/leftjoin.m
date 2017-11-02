function [t_out, Ridx, idx] = leftjoin(t_left, t_right, varargin)
% [t_out, Ridx, idx] = leftjoin(t_left, t_right, varargin)
%   outerjoin for talbe to the left, maintain order of t_left
%
%   Ridx are the indexes in t_right corresponding to the rows of t_left
%

[t_out, idx, Ridx] = outerjoin(t_left, t_right, 'type', 'left', ...
     'MergeKeys', true, varargin{:});

t_out = t_out(sortidx(idx),:);
Ridx = Ridx(sortidx(idx));

assert(all(idx>0))
warnassert(all(hist(idx,1:height(t_left))<=1), ...
    'LEFTJOIN: multiple entries from t_right match a row of t_left')
