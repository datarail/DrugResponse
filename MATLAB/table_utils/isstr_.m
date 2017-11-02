function out = isstr_( x )
%ISSTR_ return true iff argument is a "single string."
%
%    Examples:
%
%    >> isstr_('xyz')
%    ans =
%         1
%
%    >> isstr_(['x' 'y' 'z'])
%    ans =
%         1
%
%    >> isstr_(['x'; 'y'; 'z'])
%    ans =
%         0
%
%    >> ischar(['x'; 'y'; 'z'])
%    ans =
%         1
%
%    >> isstr_({'x' 'y' 'z'})
%    ans =
%         0
%
    out = ischar(x) && isrow(x);
end
