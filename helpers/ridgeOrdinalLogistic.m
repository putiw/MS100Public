function [beta, logL] = ridgeOrdinalLogistic(X, y, lambda)
% RIDGEORDINALLOGISTIC Ridge-penalized ordinal logistic regression using cumulative logit
%
% This function implements ordinal logistic regression with ridge penalty for
% ordered categorical outcomes (e.g., Low=1, Medium=2, High=3). It uses the
% cumulative logit model, which models the probability of being in category k
% or higher versus being in categories below k.
%
% Inputs:
%   X: n x p matrix of predictors (n samples, p features)
%   y: n x 1 vector of ordered categories (1, 2, ..., K)
%   lambda: scalar ridge penalty parameter (controls regularization)
%
% Outputs:
%   beta: (p+1) x 1 vector of coefficients (including intercept)
%   logL: scalar log-likelihood of the fitted model
%
% Example:
%   X = randn(100,5);  % 100 samples, 5 features
%   y = [ones(30,1); 2*ones(40,1); 3*ones(30,1)];  % 3 categories
%   lambda = 0.1;
%   [beta, logL] = ridgeOrdinalLogistic(X, y, lambda);
%
% Note: This implementation uses the cumulative logit model, which assumes
% proportional odds across categories. The ridge penalty helps prevent
% overfitting when there are many predictors.

% Initialization
X = [ones(size(X,1),1), X];  % Add intercept column
[n, p] = size(X);            % n = number of samples, p = number of features + 1 (with intercept)
K = numel(unique(y));        % K = number of categories
y = y(:);                    % Ensure y is a column vector

% Initialize parameters:
% - (K-1) threshold parameters (theta) for K categories
% - p coefficients (beta) for the predictors
theta0 = linspace(-1, 1, K-1)';  % Initial thresholds evenly spaced between -1 and 1
beta0 = zeros(p,1);              % Initial coefficients all set to zero
param0 = [theta0; beta0];        % Combine parameters into single vector

% Optimization using MATLAB's fminunc
% The negative log-likelihood function (nll_ridge_ordinal) includes:
% 1. The log-likelihood of the ordinal logistic model
% 2. The ridge penalty term (lambda * sum of squared coefficients)
opts = optimoptions('fminunc',...
    'Algorithm','quasi-newton',...  % Use quasi-Newton method
    'Display','off');              % Suppress optimization output
[param,fval] = fminunc(@(param) nll_ridge_ordinal(param, X, y, K, lambda), param0, opts);

% Extract results
beta = param(K:end);    % Get the coefficients (excluding thresholds)
logL = -fval;          % Convert negative log-likelihood to log-likelihood

end

% Note: The function nll_ridge_ordinal should be defined in a separate file
% and should implement:
% 1. The cumulative logit model probabilities
% 2. The negative log-likelihood
% 3. The ridge penalty term
%
% The cumulative logit model for category k is:
% logit(P(Y >= k)) = theta_k - X*beta
% where theta_k is the threshold for category k