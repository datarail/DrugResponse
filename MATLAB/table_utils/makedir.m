function [] = makedir(path)
    [status, msg, id] = mkdir(path);
    if status == 0
        throw(MException(id, msg));
    end
end
