function output = cellfun_dispatch(fct, varargin)

assert(isa(fct, 'function_handle'))

for i=1:length(varargin)
    if istable(varargin{i})
        varargin{i} = table2cell(varargin{i});
    end
end

varargin2 = varargin;
maxlength = max( cellfun(@length,varargin) .* (cellfun(@iscell,varargin)));
for i=1:length(varargin)
    if ~iscell(varargin{i}) || length(varargin{i})~=maxlength
        assert( ~iscell(varargin{i}) || length(varargin{2})==1 )
        
        if isvector(varargin{i}) && length(varargin{i})==maxlength
            varargin2{i} = mat2cell(varargin{i});
        else
            varargin2{i} = cell(maxlength,1);
            for j=1:maxlength
                varargin2{i}{j} = varargin{i};
            end
        end
    elseif isrow(varargin{i})
        varargin2{i} = varargin{i}';
    end        
end

output = cellfun(fct, varargin2{:}, 'uniformoutput', false);
