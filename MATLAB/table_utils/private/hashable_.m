function ktbl = hashable_(ktbl)
    w = width(ktbl);
    if w > 0 && height(ktbl) > 0
        % in principle, one row would suffice to test compatibility with
        % unique; the business with tworows below is just defensive
        % programming: guarding against possible optimizations for the
        % 1-row case in the downstream code.
        tworows = ktbl([1; 1], :);
        for i = 1:w
            try
                unique(tworows.(i));
            catch e
                if ~strcmp(e.identifier, 'MATLAB:UNIQUE:InputClass')
                    rethrow(e);
                end
                ktbl.(i) = cellmap(@DataHash, ktbl.(i));
            end
        end
    end
end
