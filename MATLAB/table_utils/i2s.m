function varargout = i2s(sz, ind)
%I2S A more convenient version of IND2SUB
%
%   If the number of requested output items is > 1, I2S(SZ, IND) should
%   behave identically to IND2SUB(SZ, IND).
%
%   Otherwise,
%
%   SUBS = I2S(SZ, IND);
%
%   ... is formally equivalent to
%
%   [I1, I2, I3, ..., In] = I2S(SZ,IND);
%   SUBS = CAT(2, I1, I2, I3, ..., In);
%
%   ...where the number of output arguments in the first assignment is
%   equal NUMEL(SZ).

    varargout = cell(1, nargout);
    if nargout > 1
        [varargout{:}] = ind2sub(sz, i);
    else
        c = cell(1, numel(sz));
        [c{:}] = ind2sub(sz, i);
        varargout{1} = cat(2, c{:});
    end
end
