function [h, p] = Modified_Mann_Kendall(x)
% Modified Mann–Kendall trend test that adjusts for autocorrelation
% Based on Hamed & Rao (1998)

x = x(:);               % Ensure column vector
n = length(x);
[~, ~, S] = Mann_Kendall(x);

% Choose a reasonable number of lags
max_lag = min(24, floor(n/4));  % Use up to 25% of length or 24 lags

% Compute autocorrelation up to max_lag
acf_vals = autocorr(x, max_lag);      % Returns (max_lag + 1)-by-1
r = acf_vals(2:end);                  % Exclude lag-0

% Ensure r is a row vector for element-wise multiplication
r = r(:).';                         

% Construct lag indices and weights
lags = 1:max_lag;
weights = 1 - lags / n;

% Compute autocorrelation correction factor (scalar)
tau = 1 + 2 * sum(weights .* r);

% Compute effective sample size
n_eff = n / tau;

% Adjusted variance of S using n_eff
varS_eff = n_eff*(n_eff-1)*(2*n_eff+5)/18;

% Compute Z statistic
if S > 0
    Z = (S - 1)/sqrt(varS_eff);
elseif S < 0
    Z = (S + 1)/sqrt(varS_eff);
else
    Z = 0;
end

% Two-tailed p-value
p = 2 * (1 - normcdf(abs(Z), 0, 1));
h = p < 0.05;
end
