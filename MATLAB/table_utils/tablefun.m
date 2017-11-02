function t_out = tablefun(fun, t_in)

t_out = t_in;
for i=1:width(t_in);
    t_out.(i) = fun(t_in.(i));
end
