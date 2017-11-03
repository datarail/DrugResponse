
function [err_tot, err_dist] = err_allrates(dist_trt, dist_0, allRates, DrugRates0, Time0, Time)
% dist_* = [tot live, tot dead, G1, S, G2, M]
%

if isvector(allRates)
    [CTM,AT] = create_TM([allRates(1:4); allRates(5:8)]);
else
    [CTM,AT] = create_TM(allRates);
end

% error for the distribution and population size
options = odeset('nonnegative', 1:6, 'abstol',1e-3, 'jacobian', [CTM zeros(4); AT zeros(4)]);
[t, test_dist_trt] = ode23(@(t,y)cell_cycle_ODE_with_dead(t,y,CTM,AT), ...
    [Time0 ToRow(Time)], [dist_0(3:end)'; ones(4,1)*dist_0(2)/4], options);
test_dist_trt = test_dist_trt(ismember(t, Time),:);
test_dist_trt = [sum(test_dist_trt(:,1:4),2) sum(test_dist_trt(:,5:8),2) test_dist_trt(:,1:4)];

% normalize error by size and sum on all time points
err_dist = (dist_trt - test_dist_trt)./max(sqrt(dist_trt),min(sqrt(dist_trt(dist_trt>0))));
err_dist = sum(err_dist,1);

err_tot = sum(err_dist([1:end 2 2]).^2) ...% more emphasis on the death count
    + sum((1e3*allRates(allRates<0)).^4) + sum(allRates(:)<0) ... % contraint for positive terms
    + sum(max(0,allRates(1,1:4)-(1.2*DrugRates0(1,1:4)))./DrugRates0(1,1:4)); % contraint on cytostatic terms;
end

