function out = sortvar(tbl, varname, order)
    col = tbl.(varname);
    if isnumeric(col) && iscell(order)
        order = cell2mat(order);
    end
    order = unique(order, 'stable');
    [~, out] = ismember(col, order);
    assert(all(out));
end

% function out = sortvar(tbl, varname, order)
%     out = zeros(height(tbl), 1, 'int16');
%
%     col = tbl.(varname);
%     col = col(:).';
%
%     if isnumeric(col) && iscell(order)
%         order = cell2mat(order);
%     end
%
%     order = unique(order, 'stable');
%
%     if iscell(col)
%         cmp = @(i) strcmp(col, order{i});
%     elseif isnumeric(col)
%         cmp = @(i) col == order(i);
%     else
%         cmp = @(i) col == order{i};
%     end
%
%     for i = 1:numel(order)
%         out(cmp(i)) = i;
%     end
% end
