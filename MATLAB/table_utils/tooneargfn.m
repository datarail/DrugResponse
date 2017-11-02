function oaf = tooneargfn(fn)
    function [varargout] = oaf_(x)
        [varargout{1:nargout}] = fn(x{:});
    end
    oaf = @oaf_;
end
