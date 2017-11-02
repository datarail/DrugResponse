function tbl = make_table(data, keyvars, valvars)
  varnames = [keyvars valvars];
    if istable(data)
        tbl = data;
        tbl.Properties.VariableNames = varnames;
    else
        tbl = table(data{:}, 'VariableNames', varnames);
        for i = 1:numel(keyvars)
            c = tbl.(keyvars{i});
            if ~iscategorical(c)
                tbl.(keyvars{i}) = ordinal(c);
                %if isnumeric(c), tbl.(keyvars{i}) = ordinal(c);
                %else tbl.(keyvars{i}) = categorical(c); end
            end
        end
    end
    userdata = make_hash({{'keyvars', keyvars}, {'valvars', valvars}});
    tbl.Properties.UserData = userdata;
end
