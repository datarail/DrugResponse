function dydt = cell_cycle_ODE_with_dead(t,y,CTM,AT)

if exist('AT','var')
    TM = CTM - AT;
else
    TM = CTM;
    AT = zeros(4);
end
   %calculate rate of change for live population (1:4) and dead ones (5:8)
   dydt = [TM * y(1:4,1);
       AT * y(1:4,1)];
end
