function tbl = make_test_table( sz, varargin )
    narginchk(1, 2);
    [sz, expanded, nd] = argchk_(sz, varargin{:});
    if nd == 0
        if expanded
            tbl = mktbl_(true, {{}, {}});
        else
            tbl = mktbl_(false, {{}}, []);
        end
        return;
    end

    n = prod(sz);
    function lvls = levels(expanded)
        lvls = make_test_ndarray(sz, expanded);
        lvls = dr.unroll(lvls, false);
        lvls = cellmap(@(i) lvls(:, i), num2cell(1:size(lvls, 2)));
    end

    keyvars = varnames_(nd, 'A');
    keylevels = levels(true);

    if expanded
        vals = keylevels;
        vns = {keyvars, varnames_(nd, 'a')};
    else
        vals = levels(false);
        vns = {keyvars};
    end

    data = [keylevels vals];
    tbl = mktbl_(expanded, vns, data{:});
end

function tbl = mktbl_(expanded, vns, varargin)
    data = varargin;
    kvs = vns{1};
    if expanded
        vvs = vns{2};
    else
        vvs = {'value'};
    end
    tbl = make_table(data, kvs, vvs);
end

function [sz, expanded, nd] = argchk_(sz, varargin)
    if iscell(sz)
        sz = cell2mat(sz);
    end
    sz = reshape(sz, 1, []);
    if ~all(sz > 0)
        error('first argument contains zero dim(s)');
    end
    nd = numel(sz);
    if nd > 26
        error('number of key variables exceeds 26');
    end
    expanded = nargin == 2 && varargin{1};
end

function vnames = varnames_(n, start)
   vnames = num2cell(char(start):char(start+n-1));
end
