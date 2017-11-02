function varargout = unpack(x, varargin)
    narginchk(1, 2);
    if nargin > 1
        d = varargin{1};
        assert(size(x, d) == nargout, ...
               sprintf(['dimension %d is not compatible ' ...
                        'with number of outputs'], d));
    else
        d = find(size(x) == nargout, 1);
        assert(~isempty(d), ...
               'no input dimension is compatible with number of outputs');
    end

    if istable(x)
        if d > 2
            throw(MException('DR20:unpack:badargs', 'invalid arguments'));
        end
        if d == 1
            tmp = rowfun(@(varargin) varargin, x);
            varargout = cell(nargout, 1);
            [varargout{:}] = unpack(tmp.(1), 1);
        else
            varargout = tabledata(x);
        end
    elseif iscell(x)
        varargout = permute(x, [d setdiff(1:ndims(x), d)]);
    else
        varargout = arraymap(@(i) hslice(x, d, i), 1:nargout);
    end

end
