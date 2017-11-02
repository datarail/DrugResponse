function out = table2mat( varargin )
    [tbl, ~, vns, ~, ~] = process_args__({'ValVars'}, varargin);
    out = table2mat_(tbl, vns);
end
