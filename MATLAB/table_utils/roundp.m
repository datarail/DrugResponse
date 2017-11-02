function out = roundp(x, p)
  if p == 0
    out = round(x);
  else
    f = 10^p;
    out = round(x*f)/f;
  end
end
