function [opt_DrugRates,final_err,exitflag] = rate_optimization_endpoint(dist_trt, dist_0, ...
    est_DrugRates, DrugRates0, Time0, Time, seed)
% dist_* = [tot live, tot dead, G1, S, G2, M]

if ~exist('RandSeed','var'), seed = round(mod(now,100)*1e5); end
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));

% seed from estimate with some randomness
seed_Rates = [min(est_DrugRates(1,:),DrugRates0(1,:)) est_DrugRates(2,:)].*(.7+.6*rand(1,8));

% adding a bit of noise and avoid values of 0 (regularization)
dist_trt = dist_trt.*(.99+.02*rand(length(Time),6));
dist_trt(:,1) = sum(dist_trt(:,3:end),2);

t0 = now;

opts = optimset('maxFunEvals', 500, 'MaxIter', 250, 'tolfun', .003, ...
    'display', 'off', 'OutputFcn', @timeup);
[opt_DrugRates,final_err,exitflag] = fminsearch( ...
    @(x) err_allrates(dist_trt, dist_0, x, DrugRates0, Time0, Time), ...
    seed_Rates, opts);

opt_DrugRates = [opt_DrugRates(1:4); opt_DrugRates(5:8)];

    function stop = timeup(x, optimValues, state)
        stop = true;
        if any(x<0)
            fprintf('-')
        elseif (now-t0)>.5/(24*60)
            fprintf('t')
        else
            stop = false;
        end
    end

end