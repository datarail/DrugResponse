function [matrix, labels] = matrix_from_table(tbl, rowvar, colvar, valvar, varargin)
%MATRIX_FROM_TABLE Extract a matrix from table.
%   M = MATRIX_FROM_TABLE(TBL, ROWVAR, COLVAR, VALVAR, AGG) extracts
%   a matrix from table TBL as follows.  First, extract the subtable
%   S = TBL(:, {ROWVAR VALVAR COLVAR}).  Then the variable COLVAR is
%   pivoted (using UNSTACK), to produce a new table T with the values
%   of S.(COLVAR) as new columns, with values derived from S.(VALVAR).
%   When more than one value S.(VALVAR) value corresponds to a single
%   entry of T, the aggregation function AGG is applied to the array of
%   all these values, and the result used as the final value in T.  If
%   AGG is omitted, it defaults to @SUM.
%
%   [M, LABELS] = MATRIX_FROM_TABLE(TBL, ROWVAR, COLVAR, VALVAR, AGG)
%   extracts into M the same matrix as described for the form above,
%   and puts in LABELS a pair of 1-dimensional tables containing
%   the row and column labels for M.

    if length(varargin) == 0
        agg = @sum;
    else
        agg = varargin{1};
    end

    warning('off', 'MATLAB:codetools:ModifiedVarnames');

    t = unstack(tbl(:, {rowvar, valvar, colvar}), ...
                valvar, colvar, 'AggregationFunction', agg);

    warning('on', 'MATLAB:codetools:ModifiedVarnames');
    
    rowlabels = t(:, rowvar);
    t.(rowvar) = [];
    if nargout > 1
      collabels = cell2table(t.Properties.VariableNames.', ...
                            'VariableNames', {colvar});
      labels = {rowlabels collabels};
    end

    % t.Properties.RowNames = rowlabels;
    % t.Properties.DimensionNames = {rowvar colvar};

    matrix = cell2mat(table2cell(t));
end
