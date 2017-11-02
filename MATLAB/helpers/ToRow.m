function row = ToRow(vect)
% row = ToRow(vect)
%   convert a vector to a row-vector

assert(any(size(vect)<=1))
if iscolumn(vect);row=vect';else row=vect;end
