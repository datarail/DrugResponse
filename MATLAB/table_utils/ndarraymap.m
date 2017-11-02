function out = ndarraymap(fn, ndarray, iterateover)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if numel(ndarray) == 0
        out = [];
        return;
    end

    n = ndims(ndarray);

    if iscell(iterateover)
        iterateover = unique(cell2mat(iterateover(:)).');
        j = numel(iterateover);
    else
        j = iterateover;
        iterateover = 1:j;
    end

    if (j >= n)
        error('bad args: too many dimensions to iterate over');
    end

    if iterateover(end) ~= j
        % if here, then iterateover is not a prefix of 1:n, so need to permute
        % the dimensions of the input n-d array
        ndarray = permute(ndarray, [iterateover setdiff(1:n, iterateover)]);
        clear('n'); % n may no longer be valid after the previous line!
    end

    out = ndarraymap_(fn, ndarray, j);

end


% function out = ndarraymap_(~, ~, ~)
% ndarray = reshape(1:prod(2:9), 2:9);
% j = 5;
% f = @(nd) numel(nd) + max(nd(:));

function out = ndarraymap_(f, ndarray, j)

    dims = size(ndarray);
    %%% m = reshape(ndarray, prod(dims(1:j)), []).';
    nr = prod(dims(1:j));
    nc = prod(dims(j+1:end));

    % reshape & cast n-d array to a nr x nc cell array
    c = mat2cell(reshape(ndarray, nr, []), ones(nr, 1), nc);

    s0 = dims(j+1:end);
    assert(size(c, 1) > 0);
    s1 = size(f(reshape(c{1}, s0)));

    if prod(s1) > 1
        s2 = [dims(1:j) s1];
    else
        s2 = dims(1:j);
    end

    tmp = cellmap(wrap_(f, s0, s1), c);
    out = reshape(cat(1, tmp{:}), s2);

end


function wrapped = wrap_(f, s0, s1)
    if s1(1) == 1
    	wrapped = @(v) f(reshape(v, s0));
    else
    	s2 = [1 prod(s1)];
    	wrapped = @(v) reshape(f(reshape(v, s0)), s2);
    end
end
