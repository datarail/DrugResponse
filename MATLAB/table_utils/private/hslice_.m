function out = hslice_(ndarray, dim, ix, subsref_template)
    subsref_template{dim} = ix;
    out = ndarray(subsref_template{:});
end
