function [grouped_tbl, starts, ends, sort_perm, sort_iperm] = group_rows_(tbl, gv)

    % The preprocessing of key columns performed by the hashable_ helper
    % function produces a table (KTBL) whose column types are among the few
    % that MATLAB's UNIQUE function can handle; KTBL is used only for
    % indexing and grouping; its contents are not included in the returned
    % value.

    ktbl = hashable_(tbl(:, dr.vidxs(tbl, gv)));
    [~, ~, group_id_col] = unique(ktbl, 'stable');
    % After the assignment
    %
    % [C, IA, group_id_col] = unique(ktbl, 'stable');
    %
    % ...the following equalities hold:
    %
    % * C(group_id_col, :) == ktbl
    % *           A(IA, :) == C

    [sorted_group_id_col, sort_perm] = sort(group_id_col);
    % Now the following equality holds:
    %
    % sorted_group_id_col == group_id_col(sort_perm)

    grouped_tbl = tbl(sort_perm, :);
    [~, sort_iperm] = sort(sort_perm);

    n = numel(sort_perm);
    assert(n == height(tbl));
    assert(isequal(sort_perm(sort_iperm).', 1:n));
    assert(isequal(sort_iperm(sort_perm).', 1:n));

    grouped_tbl.group_id = sorted_group_id_col;

    [~, starts, ~] = unique(sorted_group_id_col, 'stable');
    % After the assignment
    %
    %     [D, starts, ID] = unique(sorted_group_id_col, 'stable')
    %
    % ...the following equalities hold:
    % * sorted_grouped_id_col(starts) == D
    % *                         D(ID) == sorted_group_id_col
    %
    % wanted: vectors IE, IF such that the following equalities hold:
    %
    % * grouped_tbl(IF, :) == tbl
    % *         tbl(IE, :) == grouped_tbl

    assert(issorted(starts));
    ends = [starts(2:end); height(tbl) + 1] - 1;
end
