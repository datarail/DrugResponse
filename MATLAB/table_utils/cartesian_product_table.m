function P = cartesian_product_table(factors, factornames, varargin)
%CARTESIAN_PRODUCT_TABLE compute the cartesian product of factor groups as
%    a table.
%
%    CARTESIAN_PRODUCT_TABLE(FACTORS, FACTORNAMES) creates a table from the
%    cartesian product of FACTORS, ordered lexicographically, and with
%    column names given by FACTORNAMES.
%
%    CARTESIAN_PRODUCT_TABLE(FACTORS, FACTORNAMES, COLEX) creates a table
%    as before but orders it colexicographically if COLEX is true.

    narginchk(2, 3);
    colex = nargin > 2 && varargin{1};

    if ~(iscell(factors) && iscell(factornames))
        error('first and/or second argument is not a cell array');
    end

    nd = numel(factors);
    if numel(factornames) ~= nd
        error('numbers of factors and of factornames do not match');
    end

    if nd == 0
        P = table(0);
        P(:, 1) = []; % "empty 1-by-0 table"
    else
        content = cartesian_product(factors, colex);
        P = table(content{:, :}, 'VariableNames', factornames);
    end

end
