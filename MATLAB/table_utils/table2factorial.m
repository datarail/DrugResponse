function factorial_table = table_to_factorial(varargin)
%TABLE_TO_FACTORIAL convert arbitrary table to a factorial table.
%     FT = TABLE_TO_FACTORIAL(T) converts table T to a factorial table FT.
%
%     FT = TABLE_TO_FACTORIAL(T, 'PARAM1',val1, 'PARAM2',val2, ...) allows
%     you to specify optional parameter name/value pairs to control
%     TABLE_TO_FACTORIAL's behavior.  These parameters are listed below,
%     and their documentation is identical to those of the same-named
%     parameters for the function TABLE_TO_NDARRAY.
%
%         'KeyVars'
%         'ValVars'
%         'Aggrs'

    [tbl, kns, vns, aggrs, ~, ~] = ...
        process_args__({'KeyVars' 'ValVars' 'Aggrs'}, varargin);
    factorial_table = table_to_factorial_(tbl, kns, vns, aggrs);
end
