function out = nth(ndarray, n, dims)

%NTH gets the N-th k-dimensional subarray of NDARRAY, whose k
%   dimensions are given in DIMS
%
% [ Description ]
%   - out = nth(NDARRAY, N, DIMS)
%
% [ Arguments ]
%   - NDARRAY:  ...
%   - N:        ...
%   - DIMS:     ...
%
% [ Description ]
%   - out = nth(NDARRAY, N, DIMS) ...
%     ...
%
% [ Examples ]
%   - ...
%     \{

%     \}
%
%   - ...
%     \{

%     \}

nd = ndims(ndarray);

if iscell(dims)
    dims = unique(cell2mat(dims(:)).', 'stable');
    if max(dims) > nd || min(dims) < 1
        error('dim specs are out-of-range');
    end
    j = nd - numel(dims);
    if ~isequal(dims, j+1:nd)
        ndarray = permute(ndarray, [setdiff(1:nd, dims) dims]);
        clear('nd'); % nd may no longer be valid after the previous line!
    end
else
    j = nd - dims;
    %%% dims = (j + 1):nd;
end

sz = size(ndarray);

[idx{1:j}] = ind2sub(sz(1:j), n);
out = reshape(ndarray(idx{:}, :), sz(j+1:end));
