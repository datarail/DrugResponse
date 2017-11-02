function sentinel = make_sentinel()
%MAKE_SENTINEL make a unique sentinel constant.
%
% SENTINEL = MAKE_SENTINEL() sets SENTINEL to a unique constant value
% suitable to be used as a sentinel, by testing for identity with ISEQUAL.

    sentinel = @() 'SENTINEL';
end
