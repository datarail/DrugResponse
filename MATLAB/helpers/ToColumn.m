function col = ToColumn(vect)
% col = ToColumn(vect)
%   Convert a vector to a column vector

assert(any(size(vect)==1))
if isrow(vect);col=vect';else col=vect;end
