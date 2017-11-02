function wells = ConvertRowColToWells(rows, cols)
% wells = ConvertRowColToWells(rows, cols)

wells = strcat(char(rows+64), cellfun(@(x) {num2str(x,'%02i')}, num2cell(cols)));
