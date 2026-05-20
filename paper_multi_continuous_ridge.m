function results = paper_multi_continuous_ridge(varargin)
% PAPER_MULTI_CONTINUOUS_RIDGE Ridge linear regression analysis for continuous clinical measures
%
% Usage:
%   results = paper_multi_continuous_ridge()                           % Default: T25FW
%   results = paper_multi_continuous_ridge('measure', 'x9HPTD')        % Use x9HPTD
%   results = paper_multi_continuous_ridge('model1_metric', 'MTR')      % Different imaging metric
%
% Parameters:
%   'measure'        - Clinical measure ('T25FW', 'x9HPTD', 'x9HPTND')
%   'model1_metric'  - NRF metric ('MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm')
%   'model2_metric'  - CRF metric ('MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm')
%   'demographics'   - Demographic variables ({'Age', 'Gender', 'DurationOfDisease'})

% Returns:
%   results          - Structure containing statistical results including beta coefficients, CI, R², calibration slopes

clearvars -except varargin; clc; close all;

% Set random seed for reproducibility
rng(42);

% Parse input arguments
p = inputParser;
addParameter(p, 'measure', 'T25FW', @(x) ismember(x, {'T25FW', 'x9HPTD', 'x9HPTND'}));
addParameter(p, 'model1_metric', 'Lnorm', @(x) ismember(x, {'MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm'}));
addParameter(p, 'model2_metric', 'Lnorm', @(x) ismember(x, {'MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm'}));
addParameter(p, 'demographics', {'Age', 'Gender', 'DurationOfDisease'}, @iscell);
addParameter(p, 'alpha', 0.5, @isnumeric);
addParameter(p, 'figures', 0, @(x) ismember(x, [0, 1]));
parse(p, varargin{:});

% Extract parameters
clinical_measure = p.Results.measure;
model1_metric = p.Results.model1_metric;
model2_metric = p.Results.model2_metric;
demographic_vars = p.Results.demographics;
whichAlpha = p.Results.alpha;

% Add helper functions
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));

% Load data independently for each model to avoid excluding subjects
% that are missing data in one model but not the other
[~, ~, ~, data1, ~, y_raw1, Demo1, clinical_measure, demographic_vars] = load_ridge_data(...
    'measure', clinical_measure, ...
    'model1_metric', model1_metric, ...
    'model2_metric', model1_metric, ...
    'demographics', demographic_vars, ...
    'check_predictors', 'model1');

[~, ~, ~, ~, data2, y_raw2, Demo2, ~, ~] = load_ridge_data(...
    'measure', clinical_measure, ...
    'model1_metric', model2_metric, ...
    'model2_metric', model2_metric, ...
    'demographics', demographic_vars, ...
    'check_predictors', 'model2');

% Display configuration
fprintf('\n=== Configuration ===\n');
fprintf('Clinical Measure: %s (continuous)\n', clinical_measure);
fprintf('Model 1 (NRF) Metric: %s\n', model1_metric);
fprintf('Model 2 (CRF) Metric: %s\n', model2_metric);
fprintf('Demographics: %s\n', strjoin(demographic_vars, ', '));
fprintf('====================\n\n');

% Remove any remaining NaN or infinite values in outcome for Model 1
valid_idx1 = isfinite(y_raw1) & y_raw1 > 0;  % Ensure positive for log transform
y_raw1 = y_raw1(valid_idx1);
data1 = data1(valid_idx1, :);
if ~isempty(Demo1), Demo1 = Demo1(valid_idx1, :); end

% Log transform the outcome (often improves performance for time data)
y1 = log(y_raw1);

fprintf('MODEL 1 (NRF):\n');
fprintf('Outcome transformation: log(%s)\n', clinical_measure);
fprintf('Original range: [%.2f, %.2f]\n', min(y_raw1), max(y_raw1));
fprintf('Log-transformed range: [%.2f, %.2f]\n', min(y1), max(y1));
fprintf('Sample size after exclusions: %d\n', length(y1));
fprintf('Original %s mean: %.2f (SD: %.2f)\n', clinical_measure, mean(y_raw1), std(y_raw1));
fprintf('Log-%s mean: %.2f (SD: %.2f)\n', clinical_measure, mean(y1), std(y1));
fprintf('\n');

% Remove any remaining NaN or infinite values in outcome for Model 2
valid_idx2 = isfinite(y_raw2) & y_raw2 > 0;  % Ensure positive for log transform
y_raw2 = y_raw2(valid_idx2);
data2 = data2(valid_idx2, :);
if ~isempty(Demo2), Demo2 = Demo2(valid_idx2, :); end

% Log transform the outcome (often improves performance for time data)
y2 = log(y_raw2);

fprintf('MODEL 2 (CRF):\n');
fprintf('Outcome transformation: log(%s)\n', clinical_measure);
fprintf('Original range: [%.2f, %.2f]\n', min(y_raw2), max(y_raw2));
fprintf('Log-transformed range: [%.2f, %.2f]\n', min(y2), max(y2));
fprintf('Sample size after exclusions: %d\n', length(y2));
fprintf('Original %s mean: %.2f (SD: %.2f)\n', clinical_measure, mean(y_raw2), std(y_raw2));
fprintf('Log-%s mean: %.2f (SD: %.2f)\n', clinical_measure, mean(y2), std(y2));
fprintf('\n');

% Ridge linear regression for NRF → LP1
fprintf('Fitting NRF model with ridge regression...\n');
%%
% Standardize data
[data1_scaled, ~, ~] = zscore(data1);

lambdaGrid = logspace(-6, 6, 50);

ridgeMdl1 = fitrlinear(data1_scaled, y1, ...
    'Learner','leastsquares', ...
    'Regularization','ridge', ...
    'Lambda',lambdaGrid, ...
    'KFold',10);

cvLoss = kfoldLoss(ridgeMdl1);

% Find best Lambda
[~, bestIdx] = min(cvLoss);
bestLambda = lambdaGrid(bestIdx);

% Refit final model on all data with best Lambda
ridgeMdl1 = fitrlinear(data1_scaled, y1, ...
    'Learner','leastsquares', ...
    'Regularization','ridge', ...
    'Lambda',bestLambda);

% Extract coefficients from the final model
beta1 = ridgeMdl1.Beta;
b01 = ridgeMdl1.Bias;

% Compute linear predictors
LP1 = data1_scaled * beta1 + b01;

% Ridge linear regression for CRF → LP2
fprintf('Fitting CRF model with ridge regression...\n');

[data2_scaled, ~, ~] = zscore(data2);

finalMdl2 = fitrlinear(data2_scaled, y2, ...
    'Learner','leastsquares', ...
    'Regularization','ridge', ...
    'Lambda',lambdaGrid, ...
    'KFold',10);

cvLoss = kfoldLoss(finalMdl2);

% Find best Lambda
[~, bestIdx] = min(cvLoss);
bestLambda = lambdaGrid(bestIdx);

% Refit final model on all data with best Lambda
finalMdl2 = fitrlinear(data2_scaled, y2, ...
    'Learner','leastsquares', ...
    'Regularization','ridge', ...
    'Lambda',bestLambda);


% Extract coefficients from the final model
beta2 = finalMdl2.Beta;
b02 = finalMdl2.Bias;

% Compute linear predictors
LP2 = data2_scaled * beta2 + b02;
%% Get LP11 and LP22 from fitlm with their coefficients
fprintf('\n=== FITTING COMBINED MODELS ===\n');
mdl_combined = fitlm([LP1, Demo1], y1);
AIC_LP1 = mdl_combined.ModelCriterion.AIC;
LP11 = mdl_combined.Fitted;
coef_lp1 = mdl_combined.Coefficients.Estimate(2);  % LP1 coefficient (first predictor after intercept)
se_lp1 = mdl_combined.Coefficients.SE(2);
CI_lp1 = coef_lp1 + [-1 1] * 1.96 * se_lp1;
p_lp1 = mdl_combined.Coefficients.pValue(2);

fprintf('  LP1 Beta: %.3f [%.3f, %.3f], p = %.4f\n', coef_lp1, CI_lp1(1), CI_lp1(2), p_lp1);

mdl_combined2 = fitlm([LP2, Demo2], y2);
AIC_LP2 = mdl_combined2.ModelCriterion.AIC;
LP22 = mdl_combined2.Fitted;
coef_lp2 = mdl_combined2.Coefficients.Estimate(2);  % LP2 coefficient (first predictor after intercept)
se_lp2 = mdl_combined2.Coefficients.SE(2);
CI_lp2 = coef_lp2 + [-1 1] * 1.96 * se_lp2;
p_lp2 = mdl_combined2.Coefficients.pValue(2);

fprintf('  LP2 Beta: %.3f [%.3f, %.3f], p = %.4f\n', coef_lp2, CI_lp2(1), CI_lp2(2), p_lp2);

%% Get calibration slopes and R²
% Note: R² and calibration slopes are computed on the training data
% (apparent performance). Lambda was selected via 10-fold CV above;
% the final model is refit on all data with the best lambda, which is
% standard practice for ridge regression.
fprintf('\n=== CALIBRATION SLOPES ===\n');

% Calibration slope for LP1
mdl_cal_lp1 = fitlm(LP1, y1);
cal_slope_LP1 = mdl_cal_lp1.Coefficients.Estimate(2);
R2_LP1 = mdl_cal_lp1.Rsquared.Ordinary;
fprintf('LP1 Calibration Slope: %.3f, R²: %.3f\n', cal_slope_LP1, R2_LP1);

% Calibration slope for LP2
mdl_cal_lp2 = fitlm(LP2, y2);
cal_slope_LP2 = mdl_cal_lp2.Coefficients.Estimate(2);
R2_LP2 = mdl_cal_lp2.Rsquared.Ordinary;
fprintf('LP2 Calibration Slope: %.3f, R²: %.3f\n', cal_slope_LP2, R2_LP2);

% Calibration slope for LP11
mdl_cal_lp11 = fitlm(LP11, y1);
cal_slope_LP11 = mdl_cal_lp11.Coefficients.Estimate(2);
R2_LP11 = mdl_cal_lp11.Rsquared.Ordinary;
fprintf('LP11 Calibration Slope: %.3f, R²: %.3f\n', cal_slope_LP11, R2_LP11);

% Calibration slope for LP22
mdl_cal_lp22 = fitlm(LP22, y2);
cal_slope_LP22 = mdl_cal_lp22.Coefficients.Estimate(2);
R2_LP22 = mdl_cal_lp22.Rsquared.Ordinary;
fprintf('LP22 Calibration Slope: %.3f, R²: %.3f\n', cal_slope_LP22, R2_LP22);



fprintf('\n=== MODEL COMPARISON ===\n');
fprintf('LP1+Demo AIC: %.3f\n', AIC_LP1);
fprintf('LP2+Demo AIC: %.3f\n', AIC_LP2);

%% Compute descriptive statistics for LP1 and LP2
fprintf('\n=== DESCRIPTIVE STATISTICS ===\n');

% Descriptive statistics for LP1
LP1_min = min(LP1);
LP1_q1 = quantile(LP1, 0.25);
LP1_median = median(LP1);
LP1_q3 = quantile(LP1, 0.75);
LP1_max = max(LP1);
LP1_mean = mean(LP1);
LP1_std = std(LP1);

fprintf('LP1: min=%.3f, Q1=%.3f, median=%.3f, Q3=%.3f, max=%.3f, mean=%.3f, std=%.3f\n', ...
    LP1_min, LP1_q1, LP1_median, LP1_q3, LP1_max, LP1_mean, LP1_std);

% Descriptive statistics for LP2
LP2_min = min(LP2);
LP2_q1 = quantile(LP2, 0.25);
LP2_median = median(LP2);
LP2_q3 = quantile(LP2, 0.75);
LP2_max = max(LP2);
LP2_mean = mean(LP2);
LP2_std = std(LP2);

fprintf('LP2: min=%.3f, Q1=%.3f, median=%.3f, Q3=%.3f, max=%.3f, mean=%.3f, std=%.3f\n', ...
    LP2_min, LP2_q1, LP2_median, LP2_q3, LP2_max, LP2_mean, LP2_std);

fprintf('\n=== STORING RESULTS ===\n');

% Initialize results structure
results = struct();

results.config.clinical_measure = clinical_measure;
results.config.model1_metric = model1_metric;
results.config.model2_metric = model2_metric;
results.config.demographics = demographic_vars;
results.config.alpha = whichAlpha;

% Store LP1 results (beta from combined model)
results.LP1.beta = coef_lp1;
results.LP1.CI_lower = CI_lp1(1);
results.LP1.CI_upper = CI_lp1(2);
results.LP1.p_value = p_lp1;
results.LP1.calibration_slope = cal_slope_LP1;
results.LP1.R2 = R2_LP1;
results.LP11.R2 = R2_LP11;
results.LP1.AIC = AIC_LP1;

% Store LP1 descriptive statistics
results.LP1.min = LP1_min;
results.LP1.q1 = LP1_q1;
results.LP1.median = LP1_median;
results.LP1.q3 = LP1_q3;
results.LP1.max = LP1_max;
results.LP1.mean = LP1_mean;
results.LP1.std = LP1_std;

% Store LP2 results (beta from combined model)
results.LP2.beta = coef_lp2;
results.LP2.CI_lower = CI_lp2(1);
results.LP2.CI_upper = CI_lp2(2);
results.LP2.p_value = p_lp2;
results.LP2.calibration_slope = cal_slope_LP2;
results.LP2.R2 = R2_LP2;
results.LP22.R2 = R2_LP22;
results.LP2.AIC = AIC_LP2;

% Store LP2 descriptive statistics
results.LP2.min = LP2_min;
results.LP2.q1 = LP2_q1;
results.LP2.median = LP2_median;
results.LP2.q3 = LP2_q3;
results.LP2.max = LP2_max;
results.LP2.mean = LP2_mean;
results.LP2.std = LP2_std;

% Store LP11 results
results.LP11.calibration_slope = cal_slope_LP11;
results.LP11.R2 = R2_LP11;

% Store LP22 results
results.LP22.calibration_slope = cal_slope_LP22;
results.LP22.R2 = R2_LP22;

% Store sample sizes for each model
results.n1 = length(y1);  % Sample size for Model 1 (NRF)
results.n2 = length(y2);  % Sample size for Model 2 (CRF)

fprintf('Results stored successfully.\n');

fprintf('\n=== ANALYSIS COMPLETE ===\n');

end % Main function end