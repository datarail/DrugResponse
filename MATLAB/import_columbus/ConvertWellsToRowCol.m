function [rows, cols] = ConvertWellsToRowCol(wells)
% [rows, cols] = ConvertWellsToRowCol(wells)

rows = cellcell2cellstr(regexp(wells,'^(\w)\d*$','tokens'));
rows = [rows{:}]'-64;

cols = cellfun(@str2double, cellcell2cellstr(regexp(wells,'^\w(\d*)$','tokens')));
