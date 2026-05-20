function results_all = paper_multi_continuous_ridge_all(varargin)
% PAPER_MULTI_CONTINUOUS_RIDGE_ALL Comprehensive analysis using paper_multi_continuous_ridge.m
% Creates analysis for all continuous clinical measures with a summary table
%
% Usage:
%   results_all = paper_multi_continuous_ridge_all()                           % Default: All continuous measures
%   results_all = paper_multi_continuous_ridge_all('demographics', {'Age', 'Gender'}) % Custom demographics
%
% Parameters:
%   'demographics'   - Cell array of demographic variables ({'Age', 'Gender', 'DurationOfDisease'})
%   'figures'        - Show figures (1) or not (0) (default: 0, not currently implemented)
%
% Returns:
%   results_all      - Structure containing results for all imaging metrics and clinical measures

clearvars -except varargin; clc; close all;

% Set random seed for reproducibility
rng(42);

% Alpha parameter (passed to underlying functions)
whichAlpha = 0.05;

% Parse input arguments
p = inputParser;
addParameter(p, 'demographics', {'Age', 'Gender'}, @iscell); %, 'DurationOfDisease'
addParameter(p, 'figures', 0, @(x) ismember(x, [0, 1]));
parse(p, varargin{:});

% Extract parameters
demographic_vars = p.Results.demographics;
show_figures = p.Results.figures;

% Define clinical measures
clinical_measures = {'T25FW', 'x9HPTD', 'x9HPTND'};

% Define all imaging metrics to loop through
imaging_metrics = {'Lnorm', 'T1', 'MTR', 'FA', 'MD'};
% Alternative full set: {'LN', 'LV', 'Lnorm','MTR', 'T1', 'FA', 'MD', 'AD', 'RD'}

% Add helper functions
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));

% Display configuration
fprintf('\n=== Comprehensive Analysis Configuration (CONTINUOUS VERSION) ===\n');
fprintf('Clinical Measures: %s\n', strjoin(clinical_measures, ', '));
fprintf('Imaging Metrics: %s\n', strjoin(imaging_metrics, ', '));
fprintf('Demographics: %s\n', strjoin(demographic_vars, ', '));
fprintf('Alpha: %.2f\n', whichAlpha);
fprintf('==================================================================\n\n');

% Initialize results structure
results_all = struct();
results_all.config.demographics = demographic_vars;
results_all.config.alpha = whichAlpha;
results_all.config.clinical_measures = clinical_measures;
results_all.config.imaging_metrics = imaging_metrics;

% Initialize table data for summary
table_data = {};
row_idx = 1;

% Initialize figure data storage
figure_data = struct();

% Loop through each clinical measure and imaging metric combination
for clinical_idx = 1:length(clinical_measures)
    clinical_measure = clinical_measures{clinical_idx};
    
    for metric_idx = 1:length(imaging_metrics)
        metric = imaging_metrics{metric_idx};
        
        fprintf('=== Processing %s-%s ===\n', clinical_measure, metric);
        
        % Call paper_multi_continuous_ridge for this combination
        results_metric = paper_multi_continuous_ridge('measure', clinical_measure, ...
                                              'model1_metric', metric, ...
                                              'model2_metric', metric, ...
                                              'demographics', demographic_vars, ...
                                              'alpha', whichAlpha, ...
                                              'figures', 0);
        
        % Store results
        field_name = sprintf('%s_%s', clinical_measure, metric);
        results_all.(field_name) = results_metric;
        
        % Populate table data - one row per comparison
        % Clinical-Metric column
        table_data{row_idx, 1} = sprintf('%s-%s', clinical_measure, metric);
        table_data{row_idx, 2} = sprintf('%d/%d', results_metric.n1, results_metric.n2);  % Sample sizes (Tract/Classic)

        % Tract-based column (LP1)
        table_data{row_idx, 3} = sprintf('%.3f [%.3f,%.3f]', results_metric.LP1.beta, results_metric.LP1.CI_lower, results_metric.LP1.CI_upper);  % Beta [95% CI]
        table_data{row_idx, 4} = sprintf('%.4f', results_metric.LP1.p_value);  % p-value
        table_data{row_idx, 5} = sprintf('%.3f', results_metric.LP1.calibration_slope);  % Calibration Slope
        table_data{row_idx, 6} = sprintf('%.3f', results_metric.LP1.R2);  % R² (LP1 alone)
        table_data{row_idx, 7} = sprintf('%.3f', results_metric.LP11.R2);  % R² (LP1 + Demo)
        table_data{row_idx, 8} = sprintf('%.3f', results_metric.LP1.AIC);  % AIC

        % Classic column (LP2)
        table_data{row_idx, 9} = sprintf('%.3f [%.3f,%.3f]', results_metric.LP2.beta, results_metric.LP2.CI_lower, results_metric.LP2.CI_upper);  % Beta [95% CI]
        table_data{row_idx, 10} = sprintf('%.4f', results_metric.LP2.p_value);  % p-value
        table_data{row_idx, 11} = sprintf('%.3f', results_metric.LP2.calibration_slope);  % Calibration Slope
        table_data{row_idx, 12} = sprintf('%.3f', results_metric.LP2.R2);  % R² (LP2 alone)
        table_data{row_idx, 13} = sprintf('%.3f', results_metric.LP22.R2);  % R² (LP2 + Demo)
        table_data{row_idx, 14} = sprintf('%.3f', results_metric.LP2.AIC);  % AIC

        
        row_idx = row_idx + 1;
        
        % Store data for comprehensive figure
        figure_data.(field_name) = results_metric;
        
        fprintf('Completed %s-%s\n\n', clinical_measure, metric);
    end
end

% Create and display results table
fprintf('=== SUMMARY TABLE ===\n');
results_table = table(table_data(:,1), table_data(:,2), table_data(:,3), table_data(:,4), table_data(:,5), table_data(:,6), table_data(:,7), table_data(:,8), ...
                     table_data(:,9), table_data(:,10), table_data(:,11), table_data(:,12), table_data(:,13), table_data(:,14), ...
    'VariableNames', {'ClinicalMetric', 'N', 'TractBased_Beta_CI', 'TractBased_p_value', 'TractBased_Cal_Slope', 'TractBased_R2', 'TractBased_R2_Demo', 'TractBased_AIC', ...
                      'Classic_Beta_CI', 'Classic_p_value', 'Classic_Cal_Slope', 'Classic_R2', 'Classic_R2_Demo', 'Classic_AIC'});



% Create alternative table arrangement for Sheet 2
% Each clinical-metric combination gets 2 rows: tract-based and classic
alt_table_data = {};
alt_row_idx = 1;

% Go through the original table data and reorganize
for i = 1:length(clinical_measures)
    clinical_measure = clinical_measures{i};
    for j = 1:length(imaging_metrics)
        metric = imaging_metrics{j};

        % Find the corresponding data in results_all
        field_name = sprintf('%s_%s', clinical_measure, metric);
        if isfield(results_all, field_name)
            results_metric = results_all.(field_name);

            % Tract-based row (LP1)
            alt_table_data{alt_row_idx, 1} = sprintf('%s-%s', clinical_measure, metric);
            alt_table_data{alt_row_idx, 2} = 'tract-based';
            alt_table_data{alt_row_idx, 3} = sprintf('%d', results_metric.n1);
            alt_table_data{alt_row_idx, 4} = sprintf('%.3f [%.3f,%.3f]', results_metric.LP1.beta, results_metric.LP1.CI_lower, results_metric.LP1.CI_upper);
            alt_table_data{alt_row_idx, 5} = sprintf('%.4f', results_metric.LP1.p_value);
            alt_table_data{alt_row_idx, 6} = sprintf('%.3f', results_metric.LP1.calibration_slope);
            alt_table_data{alt_row_idx, 7} = sprintf('%.3f', results_metric.LP1.R2);
            alt_table_data{alt_row_idx, 8} = sprintf('%.3f', results_metric.LP11.R2);
            alt_table_data{alt_row_idx, 9} = sprintf('%.3f', results_metric.LP1.AIC);

            alt_row_idx = alt_row_idx + 1;

            % Classic row (LP2)
            alt_table_data{alt_row_idx, 1} = sprintf('%s-%s', clinical_measure, metric);
            alt_table_data{alt_row_idx, 2} = 'classic';
            alt_table_data{alt_row_idx, 3} = sprintf('%d', results_metric.n2);
            alt_table_data{alt_row_idx, 4} = sprintf('%.3f [%.3f,%.3f]', results_metric.LP2.beta, results_metric.LP2.CI_lower, results_metric.LP2.CI_upper);
            alt_table_data{alt_row_idx, 5} = sprintf('%.4f', results_metric.LP2.p_value);
            alt_table_data{alt_row_idx, 6} = sprintf('%.3f', results_metric.LP2.calibration_slope);
            alt_table_data{alt_row_idx, 7} = sprintf('%.3f', results_metric.LP2.R2);
            alt_table_data{alt_row_idx, 8} = sprintf('%.3f', results_metric.LP22.R2);
            alt_table_data{alt_row_idx, 9} = sprintf('%.3f', results_metric.LP2.AIC);

            alt_row_idx = alt_row_idx + 1;
        end
    end
end

% Create alternative results table
alt_results_table = table(alt_table_data(:,1), alt_table_data(:,2), alt_table_data(:,3), alt_table_data(:,4), alt_table_data(:,5), alt_table_data(:,6), alt_table_data(:,7), alt_table_data(:,8), alt_table_data(:,9), ...
    'VariableNames', {'ClinicalMetric', 'Model', 'N', 'Beta_CI', 'p_value', 'Cal_Slope', 'R2', 'R2_Demo', 'AIC'});

% Create descriptive statistics table with Sheet 2 row structure
desc_stats_data = {};
desc_row_idx = 1;

for i = 1:length(clinical_measures)
    clinical_measure = clinical_measures{i};
    for j = 1:length(imaging_metrics)
        metric = imaging_metrics{j};

        % Find the corresponding data in results_all
        field_name = sprintf('%s_%s', clinical_measure, metric);
        if isfield(results_all, field_name)
            results_metric = results_all.(field_name);

            % Tract-based row (LP1)
            desc_stats_data{desc_row_idx, 1} = sprintf('%s-%s', clinical_measure, metric);
            desc_stats_data{desc_row_idx, 2} = 'tract-based';
            desc_stats_data{desc_row_idx, 3} = sprintf('%.3f', results_metric.LP1.min);
            desc_stats_data{desc_row_idx, 4} = sprintf('%.3f', results_metric.LP1.q1);
            desc_stats_data{desc_row_idx, 5} = sprintf('%.3f', results_metric.LP1.median);
            desc_stats_data{desc_row_idx, 6} = sprintf('%.3f', results_metric.LP1.q3);
            desc_stats_data{desc_row_idx, 7} = sprintf('%.3f', results_metric.LP1.max);
            desc_stats_data{desc_row_idx, 8} = sprintf('%.3f', results_metric.LP1.mean);
            desc_stats_data{desc_row_idx, 9} = sprintf('%.3f', results_metric.LP1.std);
            desc_row_idx = desc_row_idx + 1;

            % Classic row (LP2)
            desc_stats_data{desc_row_idx, 1} = sprintf('%s-%s', clinical_measure, metric);
            desc_stats_data{desc_row_idx, 2} = 'classic';
            desc_stats_data{desc_row_idx, 3} = sprintf('%.3f', results_metric.LP2.min);
            desc_stats_data{desc_row_idx, 4} = sprintf('%.3f', results_metric.LP2.q1);
            desc_stats_data{desc_row_idx, 5} = sprintf('%.3f', results_metric.LP2.median);
            desc_stats_data{desc_row_idx, 6} = sprintf('%.3f', results_metric.LP2.q3);
            desc_stats_data{desc_row_idx, 7} = sprintf('%.3f', results_metric.LP2.max);
            desc_stats_data{desc_row_idx, 8} = sprintf('%.3f', results_metric.LP2.mean);
            desc_stats_data{desc_row_idx, 9} = sprintf('%.3f', results_metric.LP2.std);
            desc_row_idx = desc_row_idx + 1;
        end
    end
end

% Create descriptive statistics table
desc_stats_table = table(desc_stats_data(:,1), desc_stats_data(:,2), desc_stats_data(:,3), desc_stats_data(:,4), ...
    desc_stats_data(:,5), desc_stats_data(:,6), desc_stats_data(:,7), desc_stats_data(:,8), desc_stats_data(:,9), ...
    'VariableNames', {'ClinicalMetric', 'Model', 'Min', 'Q1', 'Median', 'Q3', 'Max', 'Mean', 'Std'});

% Save to Excel with multiple sheets
repo_root = fileparts(mfilename('fullpath'));
excel_filename = fullfile(repo_root, 'outputs', 'paper_multi_continuous_ridge.xlsx');
writetable(results_table, excel_filename, 'Sheet', 'Original_Layout');
writetable(alt_results_table, excel_filename, 'Sheet', 'Alternative_Layout');
writetable(desc_stats_table, excel_filename, 'Sheet', 'Descriptive_Stats');
fprintf('Excel file with three sheets saved to: %s\n', excel_filename);

% Store all tables in output structure
results_all.summary_table = results_table;
results_all.alternative_table = alt_results_table;
results_all.descriptive_stats_table = desc_stats_table;


fprintf('\n=== ANALYSIS COMPLETE ===\n');

end



