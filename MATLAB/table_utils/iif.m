function out = iif(cond, v1, v2)
%IIF stands for "inline if"
%
% [ Description ]
%   - out = iif(COND, V1, V2)
%
% [ Arguments ]
%   - COND:     conditional's test (a boolean value)
%   - V1:       value to return if COND is true
%   - V2:       value to return if COND is false
%
% [ Description ]
%   - out = iif(COND, V1, V2) returns V1 if COND is true, otherwise,
%     returns V2.
%
% [ Examples ]
%   - Conditional assignments:
%     \{
%         n = iif(mod(501, 7) == 0, 3, 4)
%        =>
%             4
%     \}
%
%   - Inline conditional expressions:
%     \{
%        ones(1, iif(mod(501, 7) == 0, 3, 4))
%        =>
%             1     1     1     1
%     \}

if cond
  out = v1;
else
  out = v2;
end
