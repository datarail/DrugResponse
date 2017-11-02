function cellArray = tsv2cell(fileName)
% cellArray = tsv2cell(fileName)
  fid = fopen(fileName,'r');   %# Open the file
  cellArray = cell(100,1);     %# Preallocate a cell array (ideally slightly
                               %#   larger than is needed)
  lineIndex = 1;               %# Index of cell to place the next line in
  nextLine = fgetl(fid);       %# Read the first line from the file
  while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
    cellArray{lineIndex} = nextLine;  %# Add the line to the cell array
    lineIndex = lineIndex+1;          %# Increment the line index
    nextLine = fgetl(fid);            %# Read the next line from the file
  end
  fclose(fid);                 %# Close the file
  cellArray = cellArray(1:lineIndex-1);  %# Remove empty cells, if needed
  for iLine = 1:lineIndex-1              %# Loop over lines
      if isempty(cellArray{iLine}), continue
      end
    lineData = textscan(cellArray{iLine},'%s',...  %# Read strings
                        'Delimiter','\t');
    lineData = lineData{1};              %# Remove cell encapsulation
    if strcmp(cellArray{iLine}(end),'\t')  %# Account for when the line
      lineData{end+1} = '';                     %#   ends with a delimiter
    end
    cellArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
  end
end
