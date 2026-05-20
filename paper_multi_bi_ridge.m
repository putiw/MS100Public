function results = paper_multi_bi_ridge(varargin)
% PAPER_MULTI_BI Ridge regression analysis for clinical measures
%
% Usage:
%   results = paper_multi_bi_ridge()                           % Default: EDSS with threshold 3
%   results = paper_multi_bi_ridge('measure', 'EDSS')         % Use EDSS with default threshold 3
%   results = paper_multi_bi_ridge('measure', 'EDSS', 'threshold', 5) % Custom threshold
%   results = paper_multi_bi_ridge('measure', 'EDSS', 'model1_metric', 'MTR', 'model2_metric', 'MTR') % Different imaging metric
%
% Parameters:
%   'measure'        - Clinical measure ('EDSS', 'T25FW', 'x9HPTD', 'x9HPTND', 'MSPro')
%   'threshold'      - Threshold for binary classification (default varies by measure)
%   'model1_metric'  - NRF metric ('MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm')
%   'model2_metric'  - CRF metric ('MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm')
%   'demographics'   - Demographic variables ({'Age', 'Gender', 'DurationOfDisease'})
%   'alpha'          - Confidence level for bootstrap CI (default: 0.05 for 95% CI)
%   'figures'        - Show figures (1) or not (0) (default: 1)
%
% Returns:
%   results          - Structure containing statistical results including OR, CI, p-values, AUC, and calibration slopes

clearvars -except varargin; clc; close all;
results = struct();
% Set random seed for reproducibility
rng(42);

% Parse input arguments
p = inputParser;
addParameter(p, 'measure', 'MSPro', @(x) ismember(x, {'EDSS', 'MSPro'}));
addParameter(p, 'threshold', [], @isnumeric);
addParameter(p, 'model1_metric', 'T1', @(x) ismember(x, {'FA', 'T1', 'MTR', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm'}));
addParameter(p, 'model2_metric', 'T1', @(x) ismember(x, {'FA', 'T1', 'MTR', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm'}));
addParameter(p, 'demographics', {'Age', 'Gender', 'DurationOfDisease'}, @iscell);
addParameter(p, 'alpha', 0.05, @isnumeric);
addParameter(p, 'figures', 1, @(x) ismember(x, [0, 1]));
parse(p, varargin{:});

% Extract parameters
clinical_measure = p.Results.measure;
model1_metric = p.Results.model1_metric;
model2_metric = p.Results.model2_metric;
demographic_vars = p.Results.demographics;
whichAlpha = p.Results.alpha;
show_figures = p.Results.figures;

% Set default thresholds if not specified
if isempty(p.Results.threshold)
    switch clinical_measure
        case 'EDSS'
            threshold = 3;
        case 'MSPro'
            threshold = 1;
        otherwise
            error('Unknown clinical measure: %s', clinical_measure);
    end
else
    threshold = p.Results.threshold;
end

% Add helper functions
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));

% Load data independently for each model to avoid excluding subjects
% that are missing data in one model but not the other
[~, ~, ~, data1, ~, y_raw1, Demo1, clinical_measure, demographic_vars] = load_ridge_data(...
    'measure', clinical_measure, ...
    'model1_metric', model1_metric, ...
    'model2_metric', model1_metric, ...
    'demographics', demographic_vars, ...
    'check_predictors', 'both');

[~, ~, ~, ~, data2, y_raw2, Demo2, ~, ~] = load_ridge_data(...
    'measure', clinical_measure, ...
    'model1_metric', model2_metric, ...
    'model2_metric', model2_metric, ...
    'demographics', demographic_vars, ...
    'check_predictors', 'both');

% Display configuration
fprintf('\n=== Configuration ===\n');
fprintf('Clinical Measure: %s (threshold: %.1f)\n', clinical_measure, threshold);
fprintf('Model 1 (NRF) Metric: %s\n', model1_metric);
fprintf('Model 2 (CRF) Metric: %s\n', model2_metric);
fprintf('Demographics: %s\n', strjoin(demographic_vars, ', '));
fprintf('====================\n\n');

% Create binary outcome variables based on threshold
y1 = y_raw1 > threshold;
y2 = y_raw2 > threshold;

% Display sample information
fprintf('MODEL 1 (NRF):\n');
fprintf('Sample size after exclusions: %d\n', length(y1));
fprintf('High %s group (>%.1f): %d (%.1f%%)\n', clinical_measure, threshold, sum(y1), 100*mean(y1));
fprintf('Low %s group (<=%.1f): %d (%.1f%%)\n', clinical_measure, threshold, sum(~y1), 100*mean(~y1));
fprintf('\n');

fprintf('MODEL 2 (CRF):\n');
fprintf('Sample size after exclusions: %d\n', length(y2));
fprintf('High %s group (>%.1f): %d (%.1f%%)\n', clinical_measure, threshold, sum(y2), 100*mean(y2));
fprintf('Low %s group (<=%.1f): %d (%.1f%%)\n', clinical_measure, threshold, sum(~y2), 100*mean(~y2));
fprintf('\n');

%%

[data1_scaled, ~, ~] = zscore(data1);

lambdaGrid = logspace(-6, 6, 50);

cvMdl = fitclinear(data1_scaled, y1, ...
    'Learner','logistic', ...
    'Regularization','ridge', ...
    'Lambda',lambdaGrid, ...
    'KFold',10);

cvLoss = kfoldLoss(cvMdl);

% Find best Lambda
[~, bestIdx] = min(cvLoss);
bestLambda = lambdaGrid(bestIdx);

% Refit final model on all data with best Lambda
Mdl1 = fitclinear(data1_scaled, y1, ...
    'Learner','logistic', ...
    'Regularization','ridge', ...
    'Lambda',bestLambda);

fprintf('bestLambda for Tract: %d\n', bestLambda);

% Extract coefficients from the final model
beta1 = Mdl1.Beta;
b01 = Mdl1.Bias;

% Compute linear predictors
LP1 = data1_scaled * beta1 + b01;


[data2_scaled, ~, ~] = zscore(data2);

cvMdl = fitclinear(data2_scaled, y2, ...
    'Learner','logistic', ...
    'Regularization','ridge', ...
    'Lambda',lambdaGrid, ...
    'KFold',10);

cvLoss = kfoldLoss(cvMdl);

% Find best Lambda
[~, bestIdx] = min(cvLoss);
bestLambda = lambdaGrid(bestIdx);

% Refit final model on all data with best Lambda
Mdl2 = fitclinear(data2_scaled, y2, ...
    'Learner','logistic', ...
    'Regularization','ridge', ...
    'Lambda',bestLambda);
    
fprintf('bestLambda for Classic: %d\n', bestLambda);


% Extract coefficients from the final model
beta2 = Mdl2.Beta;
b02 = Mdl2.Bias;

% Compute linear predictors
LP2 = data2_scaled * beta2 + b02;

% Compute descriptive statistics for LP1
results.LP1.min   = min(LP1);
results.LP1.q1    = quantile(LP1, 0.25);
results.LP1.median= median(LP1);
results.LP1.q3    = quantile(LP1, 0.75);
results.LP1.max   = max(LP1);
results.LP1.mean  = mean(LP1);
results.LP1.std   = std(LP1);

% Compute descriptive statistics for LP2
results.LP2.min   = min(LP2);
results.LP2.q1    = quantile(LP2, 0.25);
results.LP2.median= median(LP2);
results.LP2.q3    = quantile(LP2, 0.75);
results.LP2.max   = max(LP2);
results.LP2.mean  = mean(LP2);
results.LP2.std   = std(LP2);

%% Get LP11 from fitglm and OR for LP1
fprintf('\n=== FITTING COMBINED MODELS ===\n');

mdl = fitglm([LP1, Demo1], y1, 'Distribution', 'binomial');
LP11 = mdl.Fitted.LinearPredictor;
results.LP1.AIC = mdl.ModelCriterion.AIC;

coef_lp1 = mdl.Coefficients.Estimate(2);  % Coefficient for LP1
se_lp1 = mdl.Coefficients.SE(2);
OR_lp1 = exp(coef_lp1);
CI_lp1 = exp(coef_lp1 + [-1 1] * 1.96 * se_lp1);
% fitglm returns NaN pValue for binomial models (MATLAB bug with Dispersion=1)
% Use manual Wald z-test instead (equivalent, standard for logistic regression)
p_lp1 = 2 * (1 - normcdf(abs(coef_lp1 / se_lp1)));


mdl2 = fitglm([LP2, Demo2], y2, 'Distribution', 'binomial');
LP22 = mdl2.Fitted.LinearPredictor;
results.LP2.AIC = mdl2.ModelCriterion.AIC;
coef_lp2 = mdl2.Coefficients.Estimate(2);  % Coefficient for LP2
se_lp2 = mdl2.Coefficients.SE(2);
OR_lp2 = exp(coef_lp2);
CI_lp2 = exp(coef_lp2 + [-1 1] * 1.96 * se_lp2);
% fitglm returns NaN pValue for binomial models (MATLAB bug with Dispersion=1)
p_lp2 = 2 * (1 - normcdf(abs(coef_lp2 / se_lp2)));

%% Get calibration slopes
% Note: calibration slopes and AUC are computed on the training data
% (apparent performance). Lambda was selected via 10-fold CV above;
% the final model is refit on all data with the best lambda, which is
% standard practice for ridge regression.
fprintf('\n=== CALIBRATION SLOPES ===\n');

% Calibration slope for LP1
mdl_cal_lp1 = fitglm(LP1, y1, 'Distribution', 'binomial');
cal_slope_LP1 = mdl_cal_lp1.Coefficients.Estimate(2);
fprintf('LP1 Calibration Slope: %.3f\n', cal_slope_LP1);

% Calibration slope for LP2
mdl_cal_lp2 = fitglm(LP2, y2, 'Distribution', 'binomial');
cal_slope_LP2 = mdl_cal_lp2.Coefficients.Estimate(2);
fprintf('LP2 Calibration Slope: %.3f\n', cal_slope_LP2);

% Calibration slope for LP11
mdl_cal_lp11 = fitglm(LP11, y1, 'Distribution', 'binomial');
cal_slope_LP11 = mdl_cal_lp11.Coefficients.Estimate(2);
fprintf('LP11 Calibration Slope: %.3f\n', cal_slope_LP11);

% Calibration slope for LP22
mdl_cal_lp22 = fitglm(LP22, y2, 'Distribution', 'binomial');
cal_slope_LP22 = mdl_cal_lp22.Coefficients.Estimate(2);
fprintf('LP22 Calibration Slope: %.3f\n', cal_slope_LP22);
%%
%% Calculate AUC values with 95% CI using bootstrap
fprintf('\n=== AUC CALCULATION ===\n');

[AUC_nrf_mean, AUC_nrf_CI, ~] = auc_ci_bootstrap(y1, LP1, 2000, 0.05);
[AUC_crf_mean, AUC_crf_CI, ~] = auc_ci_bootstrap(y2, LP2, 2000, 0.05);
[AUC_lp11_mean, AUC_lp11_CI, ~] = auc_ci_bootstrap(y1, LP11, 2000, 0.05);
[AUC_lp22_mean, AUC_lp22_CI, ~] = auc_ci_bootstrap(y2, LP22, 2000, 0.05);

fprintf('LP1 (NRF alone): AUC=%.3f (95%% CI: %.3f–%.3f)\n', AUC_nrf_mean, AUC_nrf_CI(1), AUC_nrf_CI(2));
fprintf('LP2 (CRF alone): AUC=%.3f (95%% CI: %.3f–%.3f)\n', AUC_crf_mean, AUC_crf_CI(1), AUC_crf_CI(2));
fprintf('LP11 (LP1+Demo): AUC=%.3f (95%% CI: %.3f–%.3f)\n', AUC_lp11_mean, AUC_lp11_CI(1), AUC_lp11_CI(2));
fprintf('LP22 (LP2+Demo): AUC=%.3f (95%% CI: %.3f–%.3f)\n', AUC_lp22_mean, AUC_lp22_CI(1), AUC_lp22_CI(2));

%%
prob_LP1 = 1 ./ (1 + exp(-LP1));
prob_LP2 = 1 ./ (1 + exp(-LP2));
prob_LP11 = 1 ./ (1 + exp(-LP11));
prob_LP22 = 1 ./ (1 + exp(-LP22));

results.config.clinical_measure = clinical_measure;
results.config.threshold = threshold;
results.config.model1_metric = model1_metric;
results.config.model2_metric = model2_metric;
results.config.demographics = demographic_vars;
results.config.alpha = whichAlpha;

% Store LP1 results (OR from combined model)
results.LP1.OR = OR_lp1;
results.LP1.CI_lower = CI_lp1(1);
results.LP1.CI_upper = CI_lp1(2);
results.LP1.p_value = p_lp1;
results.LP1.calibration_slope = cal_slope_LP1;
results.LP1.AUC_mean = AUC_nrf_mean;
results.LP1.AUC_CI_lower = AUC_nrf_CI(1);
results.LP1.AUC_CI_upper = AUC_nrf_CI(2);

% Store LP2 results (OR from combined model)
results.LP2.OR = OR_lp2;
results.LP2.CI_lower = CI_lp2(1);
results.LP2.CI_upper = CI_lp2(2);
results.LP2.p_value = p_lp2;
results.LP2.calibration_slope = cal_slope_LP2;
results.LP2.AUC_mean = AUC_crf_mean;
results.LP2.AUC_CI_lower = AUC_crf_CI(1);
results.LP2.AUC_CI_upper = AUC_crf_CI(2);

% Store LP11 results
results.LP11.calibration_slope = cal_slope_LP11;
results.LP11.AUC_mean = AUC_lp11_mean;
results.LP11.AUC_CI_lower = AUC_lp11_CI(1);
results.LP11.AUC_CI_upper = AUC_lp11_CI(2);

% Store LP22 results
results.LP22.calibration_slope = cal_slope_LP22;
results.LP22.AUC_mean = AUC_lp22_mean;
results.LP22.AUC_CI_lower = AUC_lp22_CI(1);
results.LP22.AUC_CI_upper = AUC_lp22_CI(2);

% Store ROC data for replotting
results.ROC_data.LP1.probabilities = prob_LP1;
results.ROC_data.LP1.AUC_mean = AUC_nrf_mean;
results.ROC_data.LP1.AUC_CI = AUC_nrf_CI;

results.ROC_data.LP2.probabilities = prob_LP2;
results.ROC_data.LP2.AUC_mean = AUC_crf_mean;
results.ROC_data.LP2.AUC_CI = AUC_crf_CI;

results.ROC_data.LP11.probabilities = prob_LP11;
results.ROC_data.LP11.AUC_mean = AUC_lp11_mean;
results.ROC_data.LP11.AUC_CI = AUC_lp11_CI;

results.ROC_data.LP22.probabilities = prob_LP22;
results.ROC_data.LP22.AUC_mean = AUC_lp22_mean;
results.ROC_data.LP22.AUC_CI = AUC_lp22_CI;

results.ROC_data.outcome1 = y1;
results.ROC_data.outcome2 = y2;
results.ROC_data.clinical_measure = clinical_measure;
results.ROC_data.threshold = threshold;

results.n1 = length(y1);  % Sample size for Model 1 (NRF)
results.n2 = length(y2);  % Sample size for Model 2 (CRF)

fprintf('Results stored successfully.\n');

%% Generate ROC curves if requested
if show_figures
    fprintf('\n=== GENERATING ROC CURVES ===\n');
    % Pass both y1 and y2 since LP1/LP11 use y1 and LP2/LP22 use y2
    generateROCCurves(y1, y2, LP1, LP2, LP11, LP22, clinical_measure, threshold, ...
        AUC_nrf_mean, AUC_nrf_CI, AUC_crf_mean, AUC_crf_CI, ...
        AUC_lp11_mean, AUC_lp11_CI, AUC_lp22_mean, AUC_lp22_CI);
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');

end

%% Helper Functions

function generateROCCurves(y1, y2, LP1, LP2, LP11, LP22, clinical_measure, threshold, ...
    AUC_nrf_mean, AUC_nrf_CI, AUC_crf_mean, AUC_crf_CI, ...
    AUC_lp11_mean, AUC_lp11_CI, AUC_lp22_mean, AUC_lp22_CI)

% Create figure with two subplots
figure('Position', [100, 100, 1200, 500]);

% Convert linear predictors to probabilities for ROC analysis
prob_LP1 = 1 ./ (1 + exp(-LP1));
prob_LP2 = 1 ./ (1 + exp(-LP2));

% Subplot 1: LP1 vs LP2
subplot(1, 2, 1);
hold on;

% NRF curve (LP1) - using rocmetrics for proper confidence intervals
rm_nrf = rocmetrics(y1, prob_LP1, 1, 'NumBootstraps', 1000);
h_nrf = plot(rm_nrf, 'ShowConfidenceIntervals', true, 'Marker', 'none');
h_nrf.LineWidth = 2;
h_nrf.Color = [0 0.6 0];  % green
h_nrf.LineStyle = '-';

% CRF curve (LP2) - using rocmetrics for proper confidence intervals
rm_crf = rocmetrics(y2, prob_LP2, 1, 'NumBootstraps', 1000);
h_crf = plot(rm_crf, 'ShowConfidenceIntervals', true, 'Marker', 'none');
h_crf.LineWidth = 2;
h_crf.Color = [0.8 0 0];  % red
h_crf.LineStyle = '-';

% Chance diagonal
plot([0 1],[0 1],'k--','LineWidth',1.5);

% Legend
legend([h_nrf, h_crf], {sprintf('LP1 (AUC = %.3f, %.3f–%.3f)', AUC_nrf_mean, AUC_nrf_CI(1), AUC_nrf_CI(2)), ...
    sprintf('LP2 (AUC = %.3f, %.3f–%.3f)', AUC_crf_mean, AUC_crf_CI(1), AUC_crf_CI(2))}, ...
    'Location','southeast','Box','off');

% Formatting
xlabel('False Positive Rate','FontSize',12);
ylabel('True Positive Rate','FontSize',12);
title(sprintf('%s ≥ %.1f: LP1 vs LP2', clinical_measure, threshold),'FontSize',14,'FontWeight','bold');
axis square; xlim([0 1]); ylim([0 1]);
set(gca,'FontSize',11,'LineWidth',1,'TickDir','out','Box','off');
grid on; grid minor;

% Subplot 2: LP11 vs LP22
subplot(1, 2, 2);
hold on;

% Convert linear predictors to probabilities for ROC analysis
prob_LP11 = 1 ./ (1 + exp(-LP11));
prob_LP22 = 1 ./ (1 + exp(-LP22));

% LP11 curve (LP1 + Demographics) - using rocmetrics for proper confidence intervals
rm_lp11 = rocmetrics(y1, prob_LP11, 1, 'NumBootstraps', 1000);
h_lp11 = plot(rm_lp11, 'ShowConfidenceIntervals', true, 'Marker', 'none');
h_lp11.LineWidth = 2;
h_lp11.Color = [0 0.4 0];  % darker green
h_lp11.LineStyle = '-';

% LP22 curve (LP2 + Demographics) - using rocmetrics for proper confidence intervals
rm_lp22 = rocmetrics(y2, prob_LP22, 1, 'NumBootstraps', 1000);
h_lp22 = plot(rm_lp22, 'ShowConfidenceIntervals', true, 'Marker', 'none');
h_lp22.LineWidth = 2;
h_lp22.Color = [0.6 0 0];  % darker red
h_lp22.LineStyle = '-';

% Chance diagonal
plot([0 1],[0 1],'k--','LineWidth',1.5);

% Legend
legend([h_lp11, h_lp22], {sprintf('LP11 (AUC = %.3f, %.3f–%.3f)', AUC_lp11_mean, AUC_lp11_CI(1), AUC_lp11_CI(2)), ...
    sprintf('LP22 (AUC = %.3f, %.3f–%.3f)', AUC_lp22_mean, AUC_lp22_CI(1), AUC_lp22_CI(2))}, ...
    'Location','southeast','Box','off');

% Formatting
xlabel('False Positive Rate','FontSize',12);
ylabel('True Positive Rate','FontSize',12);
title(sprintf('%s ≥ %.1f: LP11 vs LP22', clinical_measure, threshold),'FontSize',14,'FontWeight','bold');
axis square; xlim([0 1]); ylim([0 1]);
set(gca,'FontSize',11,'LineWidth',1,'TickDir','out','Box','off');
grid on; grid minor;

end

function [AUC_mean, AUC_CI, AUCs] = auc_ci_bootstrap(y, s, nBoot, alpha)
% Calculate AUC and its confidence interval using bootstrap
% y: binary 0/1 vector
% s: predicted scores/probabilities (same length as y)
% nBoot: number of bootstrap resamples (e.g., 2000)
% alpha: 1 - confidence level (e.g., 0.05 for 95% CI)

if nargin < 3 || isempty(nBoot), nBoot = 2000; end
if nargin < 4 || isempty(alpha), alpha = 0.05; end

y = y(:);
s = s(:);

if numel(y) ~= numel(s), error('y and s must be same length'); end

rng(42,'twister'); % reproducibility

% Bootstrapped AUCs using bootstrp (vectorized)
AUCs = bootstrp(nBoot, @(yy,ss) perfAUC(yy,ss), y, s);

AUC_mean = mean(AUCs);
AUC_CI = quantile(AUCs, [alpha/2, 1 - alpha/2]);

end

% ---- helper: AUC via perfcurve ----
function a = perfAUC(y, s)
[~,~,~,a] = perfcurve(y, s, 1);
end
