function Wells = CorrectWellsTo0xY(Wells)

% correct well labels from \w%i to \w%02i
well_idx = cellfun(@length,Wells(:))==2;
Wells(well_idx) = strcat(cellcell2cellstr(regexp(Wells(well_idx),'^(\w)\d$','tokens')), ...
    '0',cellcell2cellstr(regexp(Wells(well_idx),'^\w(\d)$','tokens')));
