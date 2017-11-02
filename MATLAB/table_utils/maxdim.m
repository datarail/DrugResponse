function out = maxdim(x)
% >> maxdim(zeros(2, 4, 4))
% ans =
%      2
     
    [~, out] = max(size(x));
end
