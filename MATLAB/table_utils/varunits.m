function varargout = varunits(tbl, varargin)
    narginchk(1, 2);
    varargout = cell(1, nargout);
    [varargout{:}] = table_props_(tbl, 'VariableUnits', varargin{:});
end
