function varargout = varnames(tbl, varargin)
%VARNAMES Get (or set) column names.
%   The function's first argument is a table.
%
%   B = varnames(A)
%
%   ...is shorthand for
%
%   B = varnames.Properties.VariableNames
%
%
%   With a second argument, which should be a cell array of strings
%   specifying new column names, the function returns a copy of the
%   first argument with its columns renamed accordingly.  Hence,
%
%   A = varnames(A, COLNAMES)
%
%   ...simulates "renaming the columns of A".  The length of
%   COLNAMES must equal the width of A, and its contents must be
%   distinct strings suitable as column names.
%
%   If no output arguments are requested, only the
%   one-input-argument form is valid, and results in the
%   pretty-printing of the column names to the terminal.

    narginchk(1, 2);
    varargout = cell(1, nargout);
    [varargout{:}] = table_props_(tbl, 'VariableNames', varargin{:});
end
