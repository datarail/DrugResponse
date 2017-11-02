function b = logeq(x,y,precision)
% b = logeq(x,y,precision)
%   equality in the log10 domain with defined precision (default 1e-4)

if nargin<3
    precision = 1e-4;
end

b = abs(log10(x)-log10(y))<precision;

    
