function str_out = cellstr2str(cell_in, spacer)
% str_out = cellstr2str(cell_in, spacer)
%   convert a cell array of string to a single string with given spacer
%   (default = ' ')
%

if ~exist('spacer','var')
    spacer = ' ';
end

if isempty(cell_in)
    str_out = '';
    return
end

if size(cell_in,1)>1
    str_out = cell(size(cell_in,1),1);
    for i=1:size(cell_in,1)
        str_out{i} = cellstr2str(cell_in(i,:), spacer);
    end
    return
end

str_out = cell_in{1};
for i=2:length(cell_in)
    str_out = [str_out spacer cell_in{i}];
end
