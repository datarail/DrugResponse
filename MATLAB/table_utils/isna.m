function out = isna(x)
  out = isequal(x, NaN) || isequal(x, Inf) || isequal(x, '');
end
