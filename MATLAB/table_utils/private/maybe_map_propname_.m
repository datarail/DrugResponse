function prop = maybe_map_propname_(tbl, prop)
    
    map = containers.Map(fields(table().Properties), ...
                         fields(tbl.Properties));

    try
        prop = map(prop);
    catch e
        if ~strcmp(e.identifier, 'MATLAB:Containers:Map:NoKey')
            rethrow(e);
        end
    end

end