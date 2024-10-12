function [gpd_dist] = MLE_gpd_limit(extremes_vector,thres)
% create the function of the gpd pdf and cdf
gpd_pdf = @(x, k, sigma, thres) pdf("Generalized Pareto",x,-exp(k),sigma,thres);
gpd_cdf = @(x, k, sigma, thres) cdf("Generalized Pareto",x,-exp(k),sigma,thres);

% create the Max Likelihood Estimator for the function
lb = [-0.5 0 thres]; % lower bound of the dist parameter
ub = [-0.0001 Inf thres]; % upper bound of the dist parameter

%initial distribution parameters from GPD
options = statset('FunValCheck','off',MaxFunEvals=1e5,MaxIter=1e20,TolX=1e-2);

% dist_init = fitdist(extremes_vector,"GeneralizedPareto",Theta=thres);
% k_init = dist_init.k;
% sigma_init = dist_init.sigma;
% initial_params = [k_init, sigma_init, thres];

[initial_params,~] = mle(extremes_vector,pdf=gpd_pdf,cdf=gpd_cdf, start=[-0.25 std(extremes_vector) thres], Options=options);

k_init = initial_params(1);
sigma_init = initial_params(2);
initial_params(3) = thres;
 


if k_init  >= 0  % if k is more than zero, we have to limit the upper bound to ub(1) 
    sprintf('fitted shape function > 0, forced to change to upper bound %0.03f', ub(1))
    initial_params(1) = ub(1);
    ex = 1;
elseif k_init <-1/2 % if k is less than -1/2 we have to limt the lower bound of k to lb (1)
    sprintf(' initial k = %0.03f, fitted shape function < -1/2, forced to change to lower bound %0.03f', k_init, lb(1))
    initial_params(1) = lb(1);
    ex = 1;
elseif k_init<0 && k_init > ub(1) % if the -0.01 <= k < 0 where the k is small but still less than our lower bound
    initial_params(1) = ub(1);
    ex = 1;
else
    ex = 1;
end


if ex==1 % if we need to redo the MLE fitting because of k outside our target
    [param_res_est,~] = mle(extremes_vector,pdf=gpd_pdf,cdf=gpd_cdf, start=initial_params, ...
        LowerBound=lb,UpperBound=ub,Options=options);
    param_res_est(3) = thres;
    gpd_dist = makedist("GeneralizedPareto",...
        k = param_res_est(1),...
        sigma=param_res_est(2),...
        theta=param_res_est(3));

elseif ex ==0 % if we do not need to redo the MLE fitting because k is within our target
    gpd_dist = makedist("GeneralizedPareto",...
        k = k_init,...
        sigma=sigma_init,...
        theta=thres);
end


end