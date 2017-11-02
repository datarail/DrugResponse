function str = tostr_(x)
    global THEPATH;
    c = cat(2, num2cell((1:30).'), head(THEPATH, 30));
    disp(c);
%     assert(isscalar(x));
%     if isstr_(x)
%         str = x;
%     elseif iscell(x) || iscategorical(x)
%         str = cellstr(x);
%     elseif istable(x)
%         n = height(x);
%         d = 1;
%     elseif numel(x) > 0
%         sz = size(x);
%         d = find(sz > 1, 1);
%         if isempty(d)
%             d = 1;
%         end
%         n = sz(d);
%     else
%         n = 0;
%         d = 1;
%     end
%
%     for i = base:base+n-1
%         fprintf('%d\t%s\n', i, x{base+i});
%     end
end
