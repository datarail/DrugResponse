function t_data = RotatePlate180(t_data)
% t_data = RotatePlate180(t_data)

noColRow = ~all(isvariable(t_data, {'Column' 'Row'}));
if noColRow
    [Row,Column] = ConvertWellsToRowCol(cellstr(t_data.Well));
    t_data = [t_data table(Row,Column)];
end

maxrow = 2^ceil(log2(max(t_data.Row)));
maxcol = 3*2^ceil(log2(max(t_data.Column)/3));

t_data.Row = maxrow - t_data.Row +1;
t_data.Column = maxcol - t_data.Column +1;

t_data.Well = categorical(ConvertRowColToWells(t_data.Row, t_data.Column));
if noColRow
    t_data.Row = [];
    t_data.Column = [];
end
