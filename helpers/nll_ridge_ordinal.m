function nll = nll_ridge_ordinal(param, X, y, K, lambda)
% NLL_RIDGE_ORDINAL Negative log-likelihood for ridge ordinal logistic regression
%
% This function computes the negative log-likelihood of the ordinal logistic
% regression model with a ridge penalty. It implements the cumulative logit
% model for ordered categorical outcomes.
%
% Inputs:
%   param: Combined vector of [thresholds; coefficients]
%   X: n x p matrix of predictors (including intercept)
%   y: n x 1 vector of ordered categories (1, 2, ..., K)
%   K: Number of categories
%   lambda: Ridge penalty parameter
%
% Output:
%   nll: Negative log-likelihood value

% Extract parameters
theta = param(1:(K-1));     % Threshold parameters
beta = param(K:end);        % Regression coefficients

% Compute linear predictor
eta = X * beta;             % Linear predictor (n x 1)

% Initialize log-likelihood
nll = 0;

% Compute probabilities and log-likelihood for each category
for k = 1:K
    if k == 1
        % For first category: P(Y = 1) = 1 - P(Y >= 2)
        p = 1 - 1./(1 + exp(theta(1) - eta));
    elseif k == K
        % For last category: P(Y = K) = P(Y >= K)
        p = 1./(1 + exp(theta(K-1) - eta));
    else
        % For middle categories: P(Y = k) = P(Y >= k) - P(Y >= k+1)
        p = 1./(1 + exp(theta(k-1) - eta)) - 1./(1 + exp(theta(k) - eta));
    end
    
    % Add to log-likelihood for samples in this category
    idx = (y == k);
    nll = nll - sum(log(p(idx) + eps));  % Add small constant (eps) for numerical stability
end

% Add ridge penalty (excluding intercept)
ridge_penalty = lambda * sum(beta(2:end).^2);
nll = nll + ridge_penalty;

end