function table2tsv(a,filename, PrintHeaders)
% table2tsv(a, filename, PrintHeaders)
%   default: printing column headers (PrintHeaders = TRUE)

if ~exist('PrintHeaders','var')
    PrintHeaders = true;
end

cell2tsv(filename, table2cellstr(a, PrintHeaders) );
