function to = keyed_matrix_copy(from, from_labels, to, to_labels)
%KEYED_MATRIX_COPY Copy data from one key-labeled square matrix to another.
%
%   In the expression
%
%   KEYED_MATRIX_COPY(FROM, FROM_LABELS, TO, TO_LABELS)
%
%   ... FROM and TO are a square matrices whose rows (and columns) correspond to
%   the labels in FROM_LABELS and TO_LABELS, respectively.
%
%   The function returns a matrix equal to TO, after updating it with
%   corresponding values from FROM.
%
%   EXAMPLE
%
%   >> FROM = repmat(1:5.', 5, 1) .* repmat((10 .^ (0:4).'), 1, 5)
%   FROM =
%              1           2           3           4           5
%             10          20          30          40          50
%            100         200         300         400         500
%           1000        2000        3000        4000        5000
%          10000       20000       30000       40000       50000
%
%   >> FROM_LABELS = strsplit('I G E C A')
%   FROM_LABELS =
%       'I'    'G'    'E'    'C'    'A'
%
%   >> TO = zeros(7)
%   TO =
%        0     0     0     0     0     0     0
%        0     0     0     0     0     0     0
%        0     0     0     0     0     0     0
%        0     0     0     0     0     0     0
%        0     0     0     0     0     0     0
%        0     0     0     0     0     0     0
%        0     0     0     0     0     0     0
%
%   >> TO_LABELS = cellstr(('B':char('B'+6)).').'
%   TO_LABELS =
%       'B'    'C'    'D'    'E'    'F'    'G'    'H'
%
%   >> keyed_matrix_copy(FROM, FROM_LABELS, TO, TO_LABELS)
%   ans =
%        0     0     0     0     0     0     0
%        0  4000     0  3000     0  2000     0
%        0     0     0     0     0     0     0
%        0   400     0   300     0   200     0
%        0     0     0     0     0     0     0
%        0    40     0    30     0    20     0
%        0     0     0     0     0     0     0
%

    [ll, ii] = copy_idxs_(from_labels, to_labels);
    to(ii, ii) = from(ll, ll);
end

function [ll, ii] = copy_idxs_(from_labels, to_labels)
    [ll, i] = ismember(from_labels, to_labels);
    ii = i(i > 0);
end
