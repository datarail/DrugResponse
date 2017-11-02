function idx = eqtable(A,B)


hA = height(A);
hB = height(B);

assert( any([hA hB]==1) || hA==hB, ...
    'A and B should have the same size or one should be of height 1')

vNames = intersect(varnames(A), varnames(B));
if isempty(vNames)
    error('Tables A and B should have overlapping column names')
end

if hB==1 % if necessary, swap the variables such that A is the variable with one row
    hB = hA; hA = 1;
    tempA = B; B = A; A = tempA;
end

idx = true(hB,1);
% loop through the columns
for i=1:length(vNames)
    % extraction and broadcasting of the variables if a single row
    if hA==1; tempA = repmat(A.(vNames{i}), hB, 1);
    else tempA = A.(vNames{i});
    end
    tempB = B.(vNames{i});

    if iscategorical(tempA) || iscategorical(tempB) || ...
            (isnumeric(tempA) && isnumeric(tempB)) || ...
            (islogical(tempA) && islogical(tempB))
        % numeric case (much faster)
        idx = idx & (tempA==tempB);
    else
        % general case
        idx = arrayfun(@isequal, tempA, tempB) & idx;
    end
end
