function pos = strfindcell(strcell,str)
% pos = strfindcell(strcell,str)
%   find the position of a string str in a cell of strings
%   strcell.
%   return the vector pos of the same length as strcell
%   containing the pos(i) of str in strcell{i} or 0 if the
%   str is not contained in strcell
%
%

pos = zeros(size(strcell));
temp = strfind(strcell,str);
for i=1:length(temp)
    if ~isempty(temp{i})
        pos(i) = temp{i}(1);
    else
        pos(i) = 0;
    end
end
