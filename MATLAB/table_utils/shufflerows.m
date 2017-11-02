function out = shufflerows(tbl)
    out = tbl(randperm(height(tbl)), :);
end
