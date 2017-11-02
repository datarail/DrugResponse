function idx = index(item, seq, varargin)
%INDEX Return first index/indices of item in seq.
%
%    The expression
%
%    I = INDEX(A, B);
%
%    ...is almost equivalent to
%
%    [~, I] = ISMEMBER(A, B);
%
%    ...and the expression
%
%    I = INDEX(A, B, 'nostrict');
%
%    ...is almost equivalent to the sequence
%
%    [~, J] = ISMEMBER(A, B);
%    I = J(J > 0);
%    CLEAR('J');
%
%    The main differences are
%
%      - INDEX(A, B) and INDEX(A, B, 'nostrict') can be "inlined";
%      - INDEX(A, B) will raise an exception unless every entry in A is
%        present in B.
%
%    Therefore, whereas the second output argument of ISMEMBER may
%    contain zeros, the output of INDEX never does.
%

    narginchk(2, 3);
    nostrict = nargin > 2;
    if nostrict && ~isequal(varargin{1}, 'nostrict')
        throw(MException('DR20:index:badargs', 'invalid arguments'));
    end

    [~, tmp] = ismember(item, seq);
    pos = tmp > 0;
    if nostrict
        idx = tmp(pos);
    else
        if ~all(pos)
            throw(MException('DR20:index:notfound', ...
                             ['some entries in first argument where not ' ...
                              'found in second argument']));
        end
        idx = tmp;
    end
end
