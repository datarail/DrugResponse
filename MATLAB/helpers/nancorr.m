function [r,p,n] = nancorr(X,Y,type)
% [r,p,n] = nancorr(X,Y,type)
%   correlation coefficient ignoring NaN in the vectors/matrices (remove
%   all entries that are NaN in both vectors in a pairwise manner).
%
%   type is the correlation metric as defined by MATLAB
%   function can also be called as   nancorr(X, type)
%
%   r is the correlation coefficient, p is the significance and n the
%   number of samples used for calculating each pairwise correlation
%   r, p, n have a size nxn (for a single input X with n columns) or nxm
%   (for two inputs X, Y with n, resp. m, columns
%

if ~exist('type','var')
    if exist('Y','var') && isstr(Y)
        type = Y;
        clear Y
    else
        type = 'pearson';
    end
end

if exist('Y','var')
    assert(any(size(X)==size(Y)), ...
        'X and Y should be column vectors or matrices with the same number of rows')
    maxi = size(X,2);
    minj = size(X,2)+1;
    X = [X Y];

else
    maxi = size(X,2);
    minj = 2;

end

r = NaN(size(X,2));
p = r;
n = zeros(size(X,2));
for i=1:maxi
    if all(isnan(X(:,i)))
        continue
    end
    idx1 = ~isnan(X(:,i));
    n(i,i) = sum(idx1);
    r(i,i) = 1;
    if strcmp(type,'pearson')
        p(i,i)= 0;
    else
        [~,p(i,i)] = corr(X(idx1,i), X(idx1,i),'type',type);
    end
    for j=max(minj,i+1):size(X,2)
        idx = idx1 & ~isnan(X(:,j));
        if ~any(idx), continue, end
        n(i,j) = sum(idx);
        [r(i,j), p(i,j)] = corr(X(idx,i), X(idx,j),'type',type);
        n(j,i) = n(i,j);
        r(j,i) = r(i,j);
        p(j,i) = p(i,j);
    end
end

if exist('Y','var')
    r = r(1:(size(X,2)-size(Y,2)), (end-size(Y,2)+1):end);
    p = p(1:(size(X,2)-size(Y,2)), (end-size(Y,2)+1):end);
    n = n(1:(size(X,2)-size(Y,2)), (end-size(Y,2)+1):end);
end
