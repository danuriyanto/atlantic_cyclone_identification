function [distribution] = LS_GPD_limit(extremes_vector,thres)

% Empirical CDF
n = length(extremes_vector);
x_values = sort(extremes_vector);
empirical_cdf = ((1:n)-0.5)' ./ n;

% Define the GPD CDF function
gpd_cdf = @(x, params) cdf("Generalized Pareto",x,params(1),params(2),params(3));

% Objective function for least squares fitting
wgt = 1 ./ sqrt(empirical_cdf.*(1-empirical_cdf));
objective_function = @(params) sum(wgt.*(gpd_cdf(x_values, params) - empirical_cdf).^2);

% initial guess try 1
options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
initial_params = fminsearch(objective_function, [-0.01 std(extremes_vector),thres],options);
k_init = initial_params(1);
sigma_init = initial_params(2);


% Bounds for the parameters
lb = [-0.5, 0, thres]; % Lower bounds
ub = [-0.001, sigma_init, thres]; % Upper bounds

if k_init  >= 0  % if k is more than zero, we have to limit the upper bound to ub(1)
    sprintf('fitted shape function > 0, forced to change to upper bound %0.03f', ub(1))
    initial_params(1) = ub(1);
    ex = 1;
elseif k_init <-1/2 % if k is less than -1/2 we have to limit the lower bound of k to lb (1)
    sprintf(' initial k = %0.03f, fitted shape function < -1/2, forced to change to lower bound %0.03f', k_init, lb(1))
    initial_params(1) = lb(1);
    ex = 1;
elseif k_init<0 && k_init > ub(1) % if the -0.01 <= k < 0 where the k is small but still more than our lower bound
    initial_params(1) = ub(1);
    ex = 1;
else
    ex = 1;
end

if ex==1 % if we need to redo the MLE fitting because of k outside our target
    % Optimization to find the best fitting parameters using fmincon
    options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
    best_params = fmincon(objective_function, initial_params, [], [], [], [], lb, ub, [], options);
    % make the distribution object
    distribution = makedist("GeneralizedPareto",k=best_params(1),sigma=best_params(2),theta=best_params(3));
else %  if we do not need to redo the MLE fitting because k is within our target
    distribution = makedist("GeneralizedPareto",...
        k = k_init,...
        sigma=sigma_init,...
        theta=thres);
end


end
