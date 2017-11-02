function [isVar, VarIdx] = isvariable(t, VariableName)
% [isVar, VarIdx] = isvariable(t, VariableName)
%   check if a variable is in a table
%
%   VariableName is a string or a cell array
%
%   outputs:
%       - isVar (boolean arrays)
%       - VarIdx (arrays of indeces)

[isVar, VarIdx] = ismember(VariableName, t.Properties.VariableNames);
