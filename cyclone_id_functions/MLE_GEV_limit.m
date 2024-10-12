function distribution = MLE_GEV_limit(extremes_vector)
% MLE_GEV_LIMIT Fit a Generalized Extreme Value (GEV) distribution using Maximum Likelihood Estimation (MLE) with constraints.
% written by: Danu Riyanto, 2024/08/23
%
% Inputs:
%   extremes_vector - Vector containing the extreme values to be fitted.
%
% Outputs:
%   distribution - Fitted GEV distribution object.
%
% The function fits a GEV distribution with a forced negative scale parameter (sigma < 0) using MLE.
% The fitting process iterates over multiple starting points to optimize the fit based on the R-squared value.

% Set optimization options for the MLE process
opt = statset('FunValCheck','off', 'MaxFunEvals', 1e5, 'MaxIter', 1e10, 'TolX', 1e-10);

% Define the GEV PDF and CDF functions with forced negative scale parameter
gev_pdf = @(x, k, sigma, mu) pdf("Generalized Extreme Value", x, -exp(k), sigma, mu);
gev_cdf = @(x, k, sigma, mu) cdf("Generalized Extreme Value", x, -exp(k), sigma, mu);

% Set bounds for the parameters [k, sigma, mu]
lb = [-7  0  0];     % Lower bounds: k, sigma, mu
ub = [-0.7 Inf Inf]; % Upper bounds: k, sigma, mu

% Initial parameter estimation using the built-in GEV fitting function
initial_params = gevfit(extremes_vector);

% Generate a set of starting points for the parameter k within the bounds
k_start = (lb(1):0.2:ub(1))';
sigma_start = initial_params(2) * ones(size(k_start));
mu_start = initial_params(3) * ones(size(k_start));

% Create an array of initial parameter sets for the MLE process
initial_params_array = [k_start sigma_start mu_start];

% Preallocate array for R-squared values
R2 = zeros(size(k_start));

% Loop through different starting points to find the best fitting parameters
for i = 1:numel(k_start)
    % Perform MLE with constraints on the parameters
    [param_res_est, ~] = mle(extremes_vector, ...
        'pdf', gev_pdf, ...
        'cdf', gev_cdf, ...
        'start', initial_params_array(i, :), ...
        'LowerBound', lb, ...
        'UpperBound', ub, ...
        'Options', opt);
    
    % Create a distribution object for the current parameters
    iterative_distribution = makedist("GeneralizedExtremeValue", ...
                                      "k", -exp(param_res_est(1)), ...
                                      "sigma", param_res_est(2), ...
                                      "mu", param_res_est(3));
    
    % Calculate the R-squared value for the current fit
    R2(i,:) = calc_R2(extremes_vector, iterative_distribution);
end

% Identify the parameters with the highest R-squared value
[~, b] = max(R2);

% Perform a final MLE using the best initial parameters
[params, ~] = mle(extremes_vector, ...
    'pdf', gev_pdf, ...
    'cdf', gev_cdf, ...
    'start', initial_params_array(b, :), ...
    'LowerBound', lb, ...
    'UpperBound', ub, ...
    'Options', opt);

% Create the final GEV distribution object with the optimized parameters
distribution = makedist("GeneralizedExtremeValue", ...
                        "k", -exp(params(1)), ...
                        "sigma", params(2), ...
                        "mu", params(3));
end
