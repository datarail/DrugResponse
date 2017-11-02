function out = genvarname_( tbl, varargin )
%GENVARNAME_ generate a valid variable name for a table.
%
%   VARNAME = GENVARNAME_(TBL [, SEED])
%   VARNAME = GENVARNAME_(EXCL [, SEED])
%

    narginchk(1, 2);
    if nargin > 1
        seed = varargin{1};
    else
        seed = 'VAR';
    end

    if iscell(tbl)
        vs = tbl;
    else
        vs = tbl.Properties.VariableNames;
    end

    if ~verLessThan('matlab', '8.3')  % 8.3 = R2014a
        vs = matlab.lang.makeUniqueStrings([vs {seed}]);
        out = vs{end};
    else
        out = genvarname(seed, vs);
    end
end

