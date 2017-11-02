function b = isint(x)
%ISINT integer test.

    % based on http://stackoverflow.com/a/6916862/559827
    b = (x == floor(x));
end
