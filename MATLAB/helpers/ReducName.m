function Redstr = ReducName(str,RemovedCharacters,IsUp)
% Redstr = ReducName(str,RemovedCharacters,IsUp)
%
%   default RemovedCharacters = [' -/\:,' 250:3e3]
%   default IsUp = false
%

if ~exist('RemovedCharacters','var') || isempty(RemovedCharacters);
    RemovedCharacters = [' -/\:,' 250:3e3];
end

if ~exist('IsUp','var')
    IsUp = false;
end



if iscell(str)
    Redstr = cell(size(str));
    for i=1:length(str)
        Redstr{i} = ReducName(str{i},RemovedCharacters,IsUp);
    end
else
    Redstr = str;
    for i=1:length(RemovedCharacters)
        idx = find(Redstr==RemovedCharacters(i));
        Redstr(idx) = [];
    end
    if IsUp
        Redstr = upper(Redstr);
    end
end
