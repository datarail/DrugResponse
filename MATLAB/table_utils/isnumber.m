function out = isnumber(x)
  out = isnumeric(x) & ~isnan(x) & ~isinf(x);
end
