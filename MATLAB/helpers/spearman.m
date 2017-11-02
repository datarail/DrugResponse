function [r,p] = spearman(x,y)
% [r,p] = spearman(x,y)
%   overwrite of corr(x,y,'type','spearman')
%

if nargin==1
    [r,p] = corr(x,'type','spearman');
elseif nargin==2
    [r,p] = corr(x,y,'type','spearman');
end
