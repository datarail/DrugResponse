function hash = make_hash(varargin)

    narginchk(0, 1);

    if nargin == 1
        pairs = varargin{1};
        if iscell(pairs)
            get = @(i) cellmap(@(p) p{i}, pairs);
        elseif istable(pairs)
            get = @(i) table2cell(pairs(:, i));
        else
            error(['unsupported class: ' class(pairs)]);
        end
        args = arraymap(get, 1:2);
    else
        args = {};
    end

    hash = containers.Map(args{:}, 'UniformValues', false);
end

