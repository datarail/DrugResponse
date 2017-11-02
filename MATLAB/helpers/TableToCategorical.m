function t_out = TableToCategorical(t_in, varidx)
% TableToCategorical(t_in, varidx)
%   change the columns in category arrays.
%       varidx is the columns to change (either Variablenames or indices)
%       set varidx=0 for all string columns (default)
%

if ~exist('varidx','var')
    varidx=0;
end
if isnumeric(varidx)
    if varidx==0
        varidx = t_in.Properties.VariableNames;
    else
        varidx = t_in.Properties.VariableNames(varidx);
    end
end

t_out = t_in;

for ivar = ToRow(varidx)
    if ~isnumeric(t_out.(ivar{:})) && ~islogical(t_out.(ivar{:})) && ...
            ~iscategorical(t_out.(ivar{:}));
        t_out.(ivar{:}) = categorical(t_out.(ivar{:}));
    end
end
