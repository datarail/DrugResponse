function [] = tbl2tsv(path_or_handle, a, varargin)
    narginchk(2, 3);
    fdlm = char(9);
    rdlm = char(10);
    if nargin > 2; fdlm = varargin{1}; end

    function s = maybe_escape_(s)
        if isstring(s)
            pat = ['[' '"' fdlm rdlm ']'];
            if regexp(s, pat)
                s = ['"' regexprep(s, '"', '""') '"'];
            end
        elseif isscalar(s) && (isnumeric(s) || islogical(s))
            s = num2str(s);
        else
            error(['entry of an unsupported class (' class(s) ')']);
        end
    end

    if height(a) > 0
        row1 = a(1, :);
        for j = 1:width(row1)
            v = row1.(j);
            if isnumeric(v)
                a.(j) = cellstr(num2str(a.(j)));
            else
                if iscategorical(v) || ischar(v)
                    a.(j) = cellstr(a.(j));
                end
                a.(j) = cellmap(@maybe_escape_, a.(j));
            end
        end
    end

    if isa(path_or_handle, 'double') && isscalar(path_or_handle)
        fh = path_or_handle;
        doclose = false;
    elseif isstring(path_or_handle)
        fh = openfh(path_or_handle, 'Wb');
        doclose = true;
    else
        cls = class(path_or_handle);
        error(['first argument has an unexpected class (' cls ')']);
    end

    function [] = write_tbl_(a)
        for i = 1:size(a, 1)
            s = CStr2String(a(i, :), fdlm, true);
            s(end) = rdlm;
            fwrite(fh, s, 'uchar');
        end
    end

    try
        write_tbl_([cellmap(@maybe_escape_, a.Properties.VariableNames); ...
                    table2cell(a)]);
    catch exc
    end
    if doclose; fclose(fh); end
    if exist('exc', 'var'); rethrow(exc); end
end
