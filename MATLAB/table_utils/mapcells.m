function out = mapcells(cells, rule, varargin)
nvargs = length(varargin);

if nvargs > 1
  error('DR20:mapcells:TooManyArguments', ...
        'requires at most 3 arguments');
end

if isa(rule, 'containers.Map')
  if nvargs == 0
    rule = @(k) mapget(rule, k, k);
  else
    d = varargin{1};
    rule = @(k) mapget(rule, k, d);
  end
elseif nvargs == 1
  error('DR20:mapcells:InvalidThirdArgument', ...
        'accepts only 2 arguments when second one is a function handle');
end
out = cellmap(rule, cells);
