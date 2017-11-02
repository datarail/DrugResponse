function t = cell2table_withVarNames(cellstr)
% t = cell2table_withVarNames(cellstr)
%   use the first row of the cell array as variable names
%

t = cell2table(cellstr(2:end,:), 'variablenames', ...
    matlab.internal.tableUtils.makeValidName(cellstr(1,:), 'warn'));
