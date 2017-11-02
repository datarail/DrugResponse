function ith = ith_scalar(k, i)
    j = i;
    function out = ith_scalar_(k)
        if iscell(k)
            for c = reshape(k, 1, [])
                out = ith_scalar_(c{:});
                if j == 0; break; end
            end
        else
            j = j - 1;
            out = k;
        end
    end

    ith = ith_scalar_(k);

    if j > 0
        error('DR20:badindex', 'Index exceeds number of elements');
    end
end
