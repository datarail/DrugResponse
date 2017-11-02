function c = num2cellstr(x, format)
% c = num2cellstr(x, format)

if exist('format','var')
    c = cellfun(@(x) num2str(x, format), num2cell(x), 'uniformoutput',0);
else
    c = cellfun(@num2str, num2cell(x), 'uniformoutput',0);
end
