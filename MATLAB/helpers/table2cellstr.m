function cstr = table2cellstr(t_in, PrintHeaders)
% cstr = table2cellstr(t_in, PrintHeaders)
%   default: printing column headers (PrintHeaders = TRUE)

t_in = TableToString(t_in, 0);

if ~exist('PrintHeaders','var') || PrintHeaders
    cstr = [t_in.Properties.VariableNames;    table2cell(t_in)];
else
    cstr = table2cell(t_in);
end
