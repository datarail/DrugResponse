function [est_k, std_k, opt_DrugRates, min_err, exitflag] = estimate_rates_endpoint(dist_trt, dist_0, dist_ctrl, Time0, Time, Nrep)
% [est_k, std_k, opt_DrugRates, min_err, exitflag] = estimate_rates_endpoint(dist_trt, dist_0, dist_ctrl, Time0, Time, Nrep)
% dist_* = [Times] X [tot live, tot dead, G1, S, G2, M]

assert(all(size(dist_trt)==[length(Time) 6]))
assert(all(size(dist_ctrl)==[length(Time) 6]))
assert(all(size(dist_0)==[1 6]))

%% normalize to the total cells at time 0 (easier for scaling of error function)
dist_trt = dist_trt./dist_0(1);
dist_ctrl = dist_ctrl./dist_0(1);
dist_0 = dist_0/dist_0(1);

%% ensure that there is no zero values (regularize rate inversion)
RegValue = 1e-3;

dist_trt(:,3:6) = dist_trt(:,3:6)+RegValue;
dist_ctrl(:,3:6) = dist_ctrl(:,3:6)+RegValue;
dist_0(:,3:6) = dist_0(:,3:6)+RegValue;
    
%% normalize to the total cells and phase distribution (necessary)
dist_trt(:,3:6) = repmat(dist_trt(:,1),1,4).*dist_trt(:,3:6)./repmat(sum(dist_trt(:,3:6),2),1,4);
dist_ctrl(:,3:6) = repmat(dist_ctrl(:,1),1,4).*dist_ctrl(:,3:6)./repmat(sum(dist_ctrl(:,3:6),2),1,4);
dist_0(3:6) = dist_0(1)*dist_0(3:6)/sum(dist_0(3:6),2);


%% get an estiamte for the rate matrices
DrugRates0 = NaN(2,4,length(Time));
est_DrugRates = NaN(2,4,length(Time));

for i = 1:length(Time)
    [est_CTM0, est_AT0] = reverse_steady_state(dist_0, dist_ctrl(i,:), Time0, Time(i));
    DrugRates0(:,:,i) = Rates_from_CTM_AT(est_CTM0, est_AT0);
    
    [est_CTM, est_AT] = reverse_steady_state(dist_0, dist_trt(i,:), Time0, Time(i));
    est_DrugRates(:,:,i) = Rates_from_CTM_AT(est_CTM, est_AT);
end

DrugRates0 = nanmean(DrugRates0,3);
est_DrugRates = mean(est_DrugRates,3);


%% perform multiple rounds of optimization

opt_DrugRates = NaN(Nrep,8);
min_err = NaN(1,Nrep);
exitflag = NaN(1,Nrep);
parfor i = 1:Nrep
    [ks,min_err(i),exitflag(i)] = rate_optimization_endpoint(...
        dist_trt, dist_0, est_DrugRates, DrugRates0, Time0, Time, i);
    opt_DrugRates(i,:) = [ks(1,:) ks(2,:)];
end

%% find a consensus rate vector and std
logerr = max(1e-4,-log(min_err));
est_k = logerr*opt_DrugRates/sum(logerr);
std_k = std(opt_DrugRates(logerr>0,:),[],1);