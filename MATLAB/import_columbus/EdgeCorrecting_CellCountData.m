function t_corrected = EdgeCorrecting_CellCountData(t_annotated, numericfields, plotting)
% t_corrected = EdgeCorrecting_CellCountData(t_annotated, numericfields, plotting)

fprintf('Edge correction function:\n');

if ~exist('numericfields','var')
    labelfields = [{'pert_type' 'RelCellCnt' 'RelGrowth' 'DesignNumber' 'Barcode' ...
        'Untrt' 'Date' 'Row' 'Column' 'Well' 'TreatmentFile' 'Replicate', ...
        'CellLine' 'Time' 'DrugName' 'Conc'} ...
        strcat('Conc', cellfun(@(x) {num2str(x)}, num2cell(2:9)))];
    numericfields = setdiff(t_annotated.Properties.VariableNames( ...
        all(cellfun(@isnumeric, table2cell(t_annotated)))), labelfields);
    if ~isempty(numericfields)
        fprintf('\tThese numeric fields will be corrected (use ''numericfields'' to specify):\n');
        for i=1:length(numericfields)
            fprintf('\t- %s\n', numericfields{i});
        end
    end
end

if ~exist('plotting','var')
    plotting = .5; % case of plotting is major difference
end

t_plates = unique(t_annotated(:,{'Barcode','CellLine' 'Time'}));
t_corrected = t_annotated;

for ip = 1:height(t_plates)
    idx = eqtable(t_annotated,t_plates(ip,:));
    t_corrected(idx,:) = EdgeCorrection(t_annotated(idx,:), ...
        numericfields, t_plates(ip,:), plotting);
end

end


function t_out = EdgeCorrection(t_in, numericfields, t_plateinfo, plotting)

t_out = t_in;

if ~all(isvariable(t_in, {'Column' 'Row'}))
    [Row,Column] = ConvertWellsToRowCol(cellstr(t_in.Well));
    t_in = [t_in table(Row,Column)];
end



t_ctrl = t_in(t_in.pert_type=='ctl_vehicle',:);


biased = false;
if any(t_ctrl.Row==1) && any(t_ctrl.Column==1) && ...
        mod(log2(max(t_ctrl.Row)),1)==0 && mod(log2(max(t_ctrl.Column/3)),1)==0

    edge_ctrl_idx = t_ctrl.Row==1 | t_ctrl.Column==1 | t_ctrl.Row==max(t_in.Row) ...
        | t_ctrl.Column==max(t_in.Column);
    edge_idx = t_in.Row==1 | t_in.Column==1 | t_in.Row==max(t_in.Row) ...
        | t_in.Column==max(t_in.Column);

    for iFields = [{'Cellcount'} setdiff(numericfields, 'Cellcount')]
        edge_vals = t_ctrl.(iFields{:})(edge_ctrl_idx);
        center_vals = t_ctrl.(iFields{:})(~edge_ctrl_idx);
        ratio = mean(center_vals)/mean(edge_vals);
        pval = ranksum(edge_vals, center_vals);
        % correct for edges only in their is a significant difference in
        % the controls
        if pval<0.05
            fprintf('\t * Correcting for %-10s in %s (r=%.2f, p=%.3f)\n', ...
                iFields{:}(1:min(end,10)), ...
                strjoin(table2cellstr(t_plateinfo,0),'/'), ratio, pval )
            t_out.(iFields{:})(edge_idx) = t_in.(iFields{:})(edge_idx)*ratio;
            if strcmp(iFields{:}, 'Cellcount'), biased = true; end
        else
            fprintf('\t . No correc. for %-10s in %s (r=%.2f, p=%.3f)\n', ...
                iFields{:}(1:min(end,10)), ...
                strjoin(table2cellstr(t_plateinfo,0),'/'), ratio, pval )
        end
    end



else
    fprintf('\t . No correc. for  (controls not on all edges)\n', ...
        strjoin(table2cellstr(t_plateinfo,0),'/'))
end



if plotting==1 || (plotting==.5 && biased)

    plate = table2ndarray(t_in, 'keyvars', {'Row' 'Column'}, ...
        'outer', 1, 'valvars',  'Cellcount');

    get_newfigure(999, [40 100 400 300])
    get_newaxes([.1 .1 .85 .85])
    imagesc(plate(:,:,1), [0 max(max(plate(:,:,1)))])
    set(gca,'fontsize',6)
    if plotting==.5
        pause(3)
    else
        pause
    end
end
end
