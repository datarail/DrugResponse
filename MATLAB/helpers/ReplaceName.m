function Repstr = ReplaceName(str,RepalcedCharacters,NewCharacter, IsUp)
% Repstr = ReplaceName(str,RepalcedCharacters,NewCharacter, IsUp)
%
%   default Repalced = [' -/\:,' 250:3e3]
%   default NewCharacter = '_'
%   default IsUp = false
%


if ~exist('RepalcedCharacters','var') || isempty(RepalcedCharacters);
    RepalcedCharacters = [' -/\:,' 250:3e3];
end

if ~exist('NewCharacter','var') || isempty(NewCharacter);
    NewCharacter = '_';
end
if ~exist('IsUp','var')
    IsUp = false;
end

if iscell(str)
    Repstr = str;
    for i=1:length(str)
        Repstr{i} = ReplaceName(str{i},RepalcedCharacters,NewCharacter,IsUp);
    end
else
    Repstr = str;
    for i=1:length(RepalcedCharacters)
        idx = find(Repstr==RepalcedCharacters(i));
        if isempty(idx), continue, end
        Repstr(idx) = NewCharacter;
    end
    if IsUp
        Repstr = upper(Repstr);
    end

end
