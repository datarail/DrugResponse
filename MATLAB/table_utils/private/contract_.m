function out = contract_(idx)
    out = cast(str2double(strjoin(arraymap(@int2str, idx), '')), ...
               class(idx));
end
