function t_raw = DefineCellcount(t_raw, NobjField, p)

if ~isfield(p, 'Cellcount') || isempty(p.Cellcount)
    t_raw.Properties.VariableNames{NobjField{1}} = 'Cellcount';
elseif ~strcmp(p.Cellcount, 'none')
    temp = table2array(t_raw(:,NobjField));
    temp = p.Cellcount(temp);
    t_raw.Cellcount = temp;
end
if ~isfield(p, 'Deadcount') || isempty(p.Deadcount)
    t_raw.Deadcount = zeros(height(t_raw),1);
elseif ~strcmp(p.Deadcount, 'none')
    temp = table2array(t_raw(:,NobjField));
    temp = p.Deadcount(temp);
    t_raw.Deadcount = temp;
end
