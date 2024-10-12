function [distribution] = LS_GEV_limit(extremes_vector)
% LS_GEV_LIMIT Fit a Generalized Extreme Value (GEV) distribution to extreme values using least squares.
% Written by: Danu Riyanto, 2024/08/23
%
% Inputs:
%   extremes_vector - Vector containing the extreme values to be fitted.
%
% Outputs:
%   distribution - Fitted GEV distribution object.
%
% The function uses least squares fitting to match the empirical cumulative distribution function (CDF)
% of the extreme values to a GEV distribution with a force shape parameter k<0.

% Calculate the empirical CDF based on the sorted extreme values
n = length(extremes_vector);
x_values = sort(extremes_vector);
empirical_cdf = ((1:n) - 0.5)' ./ n;

% Define the GEV CDF function with parameters to be optimized
gev_cdf = @(x, params) cdf("Generalized Extreme Value", x, -exp(params(1)), params(2), params(3));

% Objective function for least squares fitting, weighted by the empirical CDF
wgt = 1 ./ sqrt(empirical_cdf .* (1 - empirical_cdf));
objective_function = @(params) sum(wgt .* (gev_cdf(x_values, params) - empirical_cdf).^2);

% Initial parameter estimation using the built-in GEV fitting function
initial_params = gevfit(extremes_vector);

% Set bounds for the parameters [k, sigma, mu]
lb = [-7  0  0]; % Lower bounds: k, sigma, mu
ub = [-0.7 Inf Inf]; % Upper bounds: k, sigma, mu

% Generate a set of starting points for the parameter k within the bounds
k_start = (lb(1):0.2:ub(1))';
sigma_start = initial_params(2) * ones(size(k_start));
mu_start = initial_params(3) * ones(size(k_start));
initial_params_array = [k_start sigma_start mu_start];

% Preallocate array for R-squared values
R2 = zeros(size(k_start));

% Optimization settings
options = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);

% Loop through different starting points to find the best fitting parameters
for i = 1:numel(k_start)
    % Perform constrained optimization to minimize the objective function
    best_params = fmincon(objective_function, initial_params_array(i, :), [], [], [], [], lb, ub, [], options);
    % Create a distribution object for the current parameters
    iterative_distribution = makedist("GeneralizedExtremeValue", ...
                                      "k", -exp(best_params(1)), ...
                                      "sigma", best_params(2), ...
                                      "mu", best_params(3));
    
    % Calculate the R-squared value for the current fit
    R2(i,:) = calc_R2(extremes_vector, iterative_distribution);
end

% Identify the parameters with the highest R-squared value
[~, b] = max(R2);

% Perform a final optimization using the best initial parameters
final_params = fmincon(objective_function, initial_params_array(b, :), [], [], [], [], lb, ub, [], options);

% Create the final GEV distribution object with the optimized parameters
distribution = makedist("GeneralizedExtremeValue", ...
                        "k", -exp(final_params(1)), ...
                        "sigma", final_params(2), ...
                        "mu", final_params(3));
end
