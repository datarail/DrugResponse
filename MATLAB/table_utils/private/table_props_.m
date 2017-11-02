function varargout = table_props_(tbl, prop, varargin)
    narginchk(2, 3);
    prop = maybe_map_propname_(tbl, prop);
    if nargin > 2
        if isstr_(varargin{1})
            tbl.Properties.(prop) = strsplit(varargin{1});
        else
            tbl.Properties.(prop) = varargin{1};
        end
        varargout = {tbl};
    else
        vals = tbl.Properties.(prop);
        if nargout == 0
            varargout = {};
            nr(vals);
        elseif nargout == 1
            varargout = {vals};
        else
            varargout = vals;
        end
    end
end
