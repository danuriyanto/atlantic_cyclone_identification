% GEV PDF and CDF check custom function
function [pdf, cdf] = gev_custom(x, k, sigma, mu)
% x: input data
% k: shape parameter
% sigma: scale parameter
% mu: location parameter

% Calculate in terms of the standard GEV distribution
y = ((x - mu) / sigma);

% Check the value of k
if k ~= 0
    % Calculate the pdf for k != 0
    pdf = exp(-((1 + k .* y) .^ (-1/k))) .* (1 + k .* y) .^ (-1/k - 1) / sigma;

    % Calculate the cdf for k != 0
    cdf = exp(-((1 + k .* y) .^ (-1/k)));
else
    % Calculate the pdf for k = 0 (Gumbel distribution)
    pdf = exp(-y - exp(-y)) / sigma;

    % Calculate the cdf for k = 0 (Gumbel distribution)
    cdf = exp(-exp(-y));
end
end
