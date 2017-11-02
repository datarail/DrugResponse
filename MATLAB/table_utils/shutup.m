function out = shutup(fn)
    warning('off','MATLAB:codetools:ModifiedVarnames')
    warning('off', 'MATLAB:table:ModifiedVarnames')
    out = fn();
    warning('on','MATLAB:codetools:ModifiedVarnames')
    warning('on', 'MATLAB:table:ModifiedVarnames')
end
