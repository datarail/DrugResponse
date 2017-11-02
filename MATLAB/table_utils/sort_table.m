function [y, varargout] = sort_table(x, varargin)
%SORT_TABLE Sort rows of table.
%
%   The motivation for implementing this function is to facilitate
%   equality-testing of tables with ISEQUAL, the typical use case being
%
%   ASSERT(ISEQUAL(SORT_TABLE(A), SORT_TABLE(B)));
%
%   (The built-in funciton SORTROWS fails when applied to tables that have
%   non-sortable columns, such as for example, columns where each
%   individual entry is a cell array.)
%
%   The two forms below
%
%   Y = sort_table(X);
%   Y = sort_table(X, COL);
%
%   ...are formally equivalent, respectively, to
%
%   Y = sortrows(cols2str(X));
%   Y = sortrows(cols2str(X), COL);
%
%   (One exception to this is when COL is such that ISEMPTY(COL) is TRUE.
%   This is explained further below.)
%
%   The forms
%
%   [Y, I] = sort_table(X);
%   [Y, I] = sort_table(X, COL);
%
%   ...are formally equivalent, respectively, to
%
%   [~, I] = sortrows(cols2str(X));
%   Y = X(I, :);
%
%   and
%
%   [~, I] = sortrows(cols2str(X), COL);
%   Y = X(I, COL);
%
%   Note that, in this case (namely, when NARGOUT == 2), the columns of the
%   returned table Y are not processed by COLS2STR, in contrast to what
%   happens when NARGOUT < 2.
%
%   The form
%
%   [Y, I] = sort_table(X, COL, 'convert');
%
%   ...is formally equivalent to
%
%   [Y, I] = sortrows(cols2str(X));
%
%   ...if ISEMPTY(COL) is TRUE, and is formally equivalent to
%
%   [Y, I] = sortrows(cols2str(X), COL);
%
%   ...if ISEMPTY(COL) is FALSE.
%
%   In forms such as
%
%   Y = sort_table(X, [], 'convert');
%   Y = sort_table(X, {}, 'convert');
%   Y = sort_table(X, '', 'convert');
%
%   ...the second and third arguments are redundant, and therefore ignored.
%   Similarly, when ISEMPTY(COL) is FALSE, the redundant third argument in
%   the form
%
%   Y = sort_table(X, COL, 'convert');
%
%   ...is ignored.  (A warning in such situations.)
%
%   Whenever ISEMPTY(COL) is TRUE for the second argument COL to
%   SORT_TABLE, the third argument must be provided, and be as shown above.
%   In other words, expressions such as
%
%   sort_table(X, [])
%   sort_table(X, {})
%   sort_table(X, '')
%
%   ...are invalid.

    % --------------------------------------------------------------------------

    narginchk(1, 3);
    nargoutchk(0, 2);

    if nargin == 2 && isempty(varargin{1}) || ...
       nargin == 3 && ~isequal(varargin{2}, 'convert')
        throw(MException('DR20:sort_table:badargs', 'invalid arguments'));
    end

    all_ = ':';
    if nargin < 2
        cols = all_;
    else
        cols = varargin{1};
        if isempty(cols)
            cols = all_;
        end
    end

    allcols = isequal(cols, all_);
    convert = nargout < 2 || nargin == 3;

    if nargout < 2 && nargin == 3
        wid = 'DR20:sort_table:redundant_args';
        if allcols
            warning(wid, 'The second and third arguments are redundant.')
        else
            warning(wid, 'The third argument is redundant.')
        end
        clear('wid');
    end

    % --------------------------------------------------------------------------

    if allcols || convert
        sortable_x = cols2str(x);
    else
        sortable_x = cols2str(x(:, cols));
        cols = all_;
    end

    out = cell(1, max(1, nargout));
    [out{:}] = sortrows(sortable_x, cols);

    assert(convert || numel(out) == 2);
    if convert
        y = out{1};
    else
        y = x(out{2}, :);
    end

    varargout = out(2:end);
end
