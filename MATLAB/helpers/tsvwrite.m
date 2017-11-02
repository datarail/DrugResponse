function tsvwrite( base_filename, M, suffix )
%TSVWRITE Write matrix M as TSV to base_filename + suffix + '.tsv'
%
%   TSVWRITE('filename.tsv', M, 'ext') writes matrix M into
%   'filename_ext.tsv'.
%
%   TSVWRITE('filename', M, 'ext') also writes matrix M into
%   'filename_ext.tsv'.
%
%   Uses unix newlines and enough precision for 64-bit doubles. Designed as
%   a drop-in replacement for xlswrite(FILE,ARRAY,SHEET), with the sheets
%   stored as separate suffixed .tsv files. Accepts mixed numeric/string
%   cell arrays.

% Ensure M is a 2-dimensional numeric or cell array.
assert(ismatrix(M), 'M must be 2-dimensional.');
assert(isnumeric(M) || iscell(M), 'M must be a numeric or cell array.');

% Append suffix to base_filename.
filename = [base_filename '_' suffix '.tsv'];

% Normalize and validate data.
if isnumeric(M)
    % If passed a numeric matrix, convert to cell array.
    M = num2cell(M);
else
    % If passed a cell array, ensure it only contains numbers and strings.
    assert(all(all(cellfun(@(v) isnumeric(v) || ischar(v), M))), ...
        'Cell array may only contain numeric and string values.');
end

fid = fopen(filename, 'w');
[num_rows, num_columns] = size(M);
tab = char(9);
cr = char(10);
dquote = '"';
for row = 1:num_rows
    for column = 1:num_columns
        value = M{row,column};
        if isnumeric(value)
            fprintf(fid, '%.17g', value);
        else
            % Perform escaping and quoting.
            value = strrep(value, '"', '""');
            if ~( ...
                    isempty(strfind(value, tab)) && ...
                    isempty(strfind(value, cr)) && ...
                    isempty(strfind(value, dquote)) ...
                    )
                value = ['"' value '"'];
            end
            fprintf(fid, '%s', value);
        end
        if column < num_columns
            fprintf(fid, '\t');
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);

end
