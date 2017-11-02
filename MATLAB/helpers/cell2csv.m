function cell2csv(filename,data,delimiter)
%  cell2csv(filename,data,delimiter)
%       default delimiter = ','


if ~exist('delimiter','var')
    delimiter = ',';
end
if delimiter==' '
    error('not working, needs to be corrected')
end

numidx = cellfun(@isnumeric,data(:));
data(numidx) = cellfun(@num2str,data(numidx),'uniformoutput',0);

logicidx = cellfun(@islogical,data(:));
data(logicidx) = cellfun(@num2str,data(logicidx),'uniformoutput',0);



file = fopen(filename,'w');
for i=1:size(data,1)
    idx = strfindcell(data(i,:),'\')>0;
    if any(idx)
        for j = find(idx)
            data{i,j} = strjoin(strcat(regexp(data{i,j},'\','split'),'\\'),'');
        end
    end

    data(i,1:end-1) = strcat(data(i,1:end-1),delimiter);
    fprintf(file,[data{i,:}]);
    if i<size(data,1)
        fprintf(file,'\n');
    end
end

pause(.1)
fclose(file);
