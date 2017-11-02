function [CTM, AT] = reverse_steady_state(dist_0, dist_end, Time)
% dist_* = [tot live, tot dead, G1, S, G2, M]
%
% AT estimated as the same death rate for all phases

div_rate = (1 + max(0,dist_end(2)-dist_0(2))/(dist_end(1)-dist_0(1)) )*log(dist_end(1)./dist_0(1))/Time;
% needs to be checked
div_rate = max(0,div_rate);
death_rate = max(0, ((dist_end(2)-dist_0(2))/(dist_end(1)-dist_0(1)) )*log(dist_end(1)./dist_0(1)))/Time;

% split the death rate among all phases

TM = zeros(4);
TM(1,4) = 2*div_rate * dist_end(1)/dist_end(6);
TM(4,4) = -div_rate * dist_end(1)/dist_end(6) - death_rate;

TM(4,3) = (div_rate-TM(4,4)) * dist_end(6)/dist_end(5);
TM(3,3) = -TM(4,3) -death_rate;

TM(3,2) = (div_rate-TM(3,3)) * dist_end(5)/dist_end(4);
TM(2,2) = -TM(3,2) -death_rate;

TM(2,1) = (div_rate-TM(2,2)) * dist_end(4)/dist_end(3);
TM(1,1) = -TM(2,1) -death_rate;

AT = eye(4) * death_rate;
CTM = TM + AT;