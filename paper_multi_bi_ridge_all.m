function results_all = paper_multi_bi_ridge_all(varargin)
% this is for the multi variate analysis for the binary cases like EDSS and
% MsPro
% Creates analysis for both EDSS and MSPro with a 4x5 subplot figure
% Uses the dev version that gets LP11/LP22 from fitglm and OR for LP1/LP2 from combined models
%
% Usage:
%   results_all = paper_multi_bi_dev_all()                           % Default: Both EDSS and MSPro
%   results_all = paper_multi_bi_dev_all('demographics', {'Age', 'Gender'}) % Custom demographics
%   results_all = paper_multi_bi_dev_all('alpha', 0.1)               % Custom alpha
%
% Parameters:
%   'demographics'   - Cell array of demographic variables ({'Age', 'Gender', 'DurationOfDisease'})
%   'alpha'          - Ridge regression alpha parameter (default: 0.5)
%   'figures'        - Show figures (1) or not (0) (default: 1)
%
% Returns:
%   results_all      - Structure containing results for all imaging metrics and clinical measures

clearvars -except varargin; clc; close all;

% Set random seed for reproducibility
rng(42);

% Parse input arguments
p = inputParser;
addParameter(p, 'demographics', {'Age', 'Gender', 'DurationOfDisease'}, @iscell); %, 'DurationOfDisease'
addParameter(p, 'alpha', 0.5, @isnumeric);
addParameter(p, 'figures', 0, @(x) ismember(x, [0, 1]));
parse(p, varargin{:});

% Extract parameters
demographic_vars = p.Results.demographics;
whichAlpha = p.Results.alpha;
show_figures = p.Results.figures;

% Define clinical measures and their thresholds
clinical_measures = {'EDSS', 'MSPro'};
thresholds = [3, 1];  % EDSS threshold 3, MSPro threshold 1

% Define all imaging metrics to loop through
imaging_metrics = {'Lnorm', 'T1', 'MTR', 'FA', 'MD'};

% Add helper functions
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));

% Display configuration
fprintf('\n=== Analysis Configuration (DEV VERSION) ===\n');
fprintf('Clinical Measures: %s\n', strjoin(clinical_measures, ', '));
fprintf('Imaging Metrics: %s\n', strjoin(imaging_metrics, ', '));
fprintf('Demographics: %s\n', strjoin(demographic_vars, ', '));
fprintf('Alpha: %.2f\n', whichAlpha);
fprintf('===========================================================\n\n');

% Initialize results structure
results_all = struct();
results_all.config.demographics = demographic_vars;
results_all.config.alpha = whichAlpha;
results_all.config.clinical_measures = clinical_measures;
results_all.config.imaging_metrics = imaging_metrics;
results_all.config.thresholds = thresholds;

% Initialize table data for summary
table_data = {};
row_idx = 1;

% Initialize figure data storage
figure_data = struct();

% Loop through each clinical measure and imaging metric combination
for clinical_idx = 1:length(clinical_measures)
    clinical_measure = clinical_measures{clinical_idx};
    threshold = thresholds(clinical_idx);
    
    for metric_idx = 1:length(imaging_metrics)
        metric = imaging_metrics{metric_idx};
        
        fprintf('=== Processing %s-%s ===\n', clinical_measure, metric);
        
        % Call paper_multi_bi_dev for this combination
        results_metric = paper_multi_bi_ridge('measure', clinical_measure, ...
                                          'threshold', threshold, ...
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

        table_data{row_idx, 3} = sprintf('%.2f [%.2f,%.2f]', results_metric.LP1.OR, results_metric.LP1.CI_lower, results_metric.LP1.CI_upper);  % OR [95% CI]
        table_data{row_idx, 4} = sprintf('%.4f', results_metric.LP1.p_value);  % p-value
        table_data{row_idx, 5} = sprintf('%.3f', results_metric.LP1.calibration_slope);  % Calibration Slope
        table_data{row_idx, 6} = sprintf('%.3f [%.3f,%.3f]', results_metric.LP1.AUC_mean, results_metric.LP1.AUC_CI_lower, results_metric.LP1.AUC_CI_upper);  % AUC [95% CI]
        table_data{row_idx, 7} = sprintf('%.3f', results_metric.LP1.AIC);  % AIC


        table_data{row_idx, 8} = sprintf('%.2f [%.2f,%.2f]', results_metric.LP2.OR, results_metric.LP2.CI_lower, results_metric.LP2.CI_upper);  % OR [95% CI]
        table_data{row_idx, 9} = sprintf('%.4f', results_metric.LP2.p_value);  % p-value
        table_data{row_idx, 10} = sprintf('%.3f', results_metric.LP2.calibration_slope);  % Calibration Slope
        table_data{row_idx, 11} = sprintf('%.3f [%.3f,%.3f]', results_metric.LP2.AUC_mean, results_metric.LP2.AUC_CI_lower, results_metric.LP2.AUC_CI_upper);  % AUC [95% CI]
        table_data{row_idx, 12} = sprintf('%.3f', results_metric.LP2.AIC);  % AIC
        
        row_idx = row_idx + 1;
        
        % Store ROC data for comprehensive figure
        figure_data.(field_name) = results_metric.ROC_data;
        
        fprintf('Completed %s-%s\n\n', clinical_measure, metric);
    end
end

% Create and display results table
fprintf('=== SUMMARY TABLE ===\n');
results_table = table(table_data(:,1), table_data(:,2), table_data(:,3), table_data(:,4), table_data(:,5), table_data(:,6), table_data(:,7), table_data(:,8), ...
                     table_data(:,9), table_data(:,10), table_data(:,11), table_data(:,12), ...
    'VariableNames', {'ClinicalMetric', 'N', 'TractBased_OR_CI', 'TractBased_p_value', 'TractBased_Cal_Slope', 'TractBased_AUC_CI', 'TractBased_AIC', ...
                      'Classic_OR_CI', 'Classic_p_value', 'Classic_Cal_Slope', 'Classic_AUC_CI', 'Classic_AIC'});



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
            alt_table_data{alt_row_idx, 4} = sprintf('%.2f [%.2f,%.2f]', results_metric.LP1.OR, results_metric.LP1.CI_lower, results_metric.LP1.CI_upper);
            alt_table_data{alt_row_idx, 5} = sprintf('%.4f', results_metric.LP1.p_value);
            alt_table_data{alt_row_idx, 6} = sprintf('%.3f', results_metric.LP1.calibration_slope);
            alt_table_data{alt_row_idx, 7} = sprintf('%.3f [%.3f,%.3f]', results_metric.LP11.AUC_mean, results_metric.LP11.AUC_CI_lower, results_metric.LP11.AUC_CI_upper);
            alt_table_data{alt_row_idx, 8} = sprintf('%.3f', results_metric.LP1.AIC);

            alt_row_idx = alt_row_idx + 1;

            % Classic row (LP2)
            alt_table_data{alt_row_idx, 1} = sprintf('%s-%s', clinical_measure, metric);
            alt_table_data{alt_row_idx, 2} = 'classic';
            alt_table_data{alt_row_idx, 3} = sprintf('%d', results_metric.n2);
            alt_table_data{alt_row_idx, 4} = sprintf('%.2f [%.2f,%.2f]', results_metric.LP2.OR, results_metric.LP2.CI_lower, results_metric.LP2.CI_upper);
            alt_table_data{alt_row_idx, 5} = sprintf('%.4f', results_metric.LP2.p_value);
            alt_table_data{alt_row_idx, 6} = sprintf('%.3f', results_metric.LP2.calibration_slope);
            alt_table_data{alt_row_idx, 7} = sprintf('%.3f [%.3f,%.3f]', results_metric.LP22.AUC_mean, results_metric.LP22.AUC_CI_lower, results_metric.LP22.AUC_CI_upper);
            alt_table_data{alt_row_idx, 8} = sprintf('%.3f', results_metric.LP2.AIC);
            alt_row_idx = alt_row_idx + 1;
        end
    end
end

% Create alternative results table
alt_results_table = table(alt_table_data(:,1), alt_table_data(:,2), alt_table_data(:,3), alt_table_data(:,4), alt_table_data(:,5), alt_table_data(:,6), alt_table_data(:,7), alt_table_data(:,8), ...
    'VariableNames', {'ClinicalMetric', 'Model', 'N', 'OR_CI', 'p_value', 'Cal_Slope', 'AUC_CI', 'AIC'});

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
excel_filename = fullfile(repo_root, 'outputs', 'paper_multi_bi_ridge.xlsx');
writetable(results_table, excel_filename, 'Sheet', 'Original_Layout');
writetable(alt_results_table, excel_filename, 'Sheet', 'Alternative_Layout');
writetable(desc_stats_table, excel_filename, 'Sheet', 'Descriptive_Stats');
fprintf('Excel file with three sheets saved to: %s\n', excel_filename);

% Store all tables in output structure
results_all.summary_table = results_table;
results_all.alternative_table = alt_results_table;
results_all.descriptive_stats_table = desc_stats_table;

% Generate  figure if requested
if show_figures
    fprintf('\n=== GENERATING FIGURE ===\n');
    generateFigure(figure_data, imaging_metrics);
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');

end

%% Helper Functions



function generateFigure(figure_data, imaging_metrics)
% Generate 4x5 subplot figure with ROC curves and AUC (95% CI) displayed on each subplot
% Row 1: LP1 and LP2 for EDSS (5 metrics)
% Row 2: LP11 and LP22 for EDSS (5 metrics) 
% Row 3: LP1 and LP2 for MSPro (5 metrics)
% Row 4: LP11 and LP22 for MSPro (5 metrics)

% Create figure
figure(1); clf;
set(gcf, 'Position', [100, 100, 2000, 1200]);

% Define colors
color_lp1 = [0 0.6 0];      % Green for LP1
color_lp2 = [0.8 0 0];      % Red for LP2
color_lp11 = [0 0.4 0];     % Darker green for LP11
color_lp22 = [0.6 0 0];     % Darker red for LP22

% Row 1: LP1 and LP2 for EDSS
for metric_idx = 1:length(imaging_metrics)
    metric = imaging_metrics{metric_idx};
    field_name = sprintf('EDSS_%s', metric);
    
    if isfield(figure_data, field_name)
        roc_data = figure_data.(field_name);
        subplot(4, 5, metric_idx);
        hold on;
        
        % Extract data
        y = roc_data.outcome1;  % Use outcome1 for LP1
        prob_LP1 = roc_data.LP1.probabilities;
        prob_LP2 = roc_data.LP2.probabilities;
        
        % Plot ROC curves with confidence intervals
        % LP1 curve
        rm_lp1 = rocmetrics(y, prob_LP1, 1, 'NumBootstraps', 1000);
        h_lp1 = plot(rm_lp1, 'ShowConfidenceIntervals', true, 'Marker', 'none');
        h_lp1.LineWidth = 2;
        h_lp1.Color = color_lp1;
        h_lp1.LineStyle = '-';
        
        % LP2 curve
        rm_lp2 = rocmetrics(y, prob_LP2, 1, 'NumBootstraps', 1000);
        h_lp2 = plot(rm_lp2, 'ShowConfidenceIntervals', true, 'Marker', 'none');
        h_lp2.LineWidth = 2;
        h_lp2.Color = color_lp2;
        h_lp2.LineStyle = '-';
        
        % Chance diagonal
        plot([0 1],[0 1],'k--','LineWidth',1.5);
        
        % Add AUC (95% CI) text on the plot
        auc_text = sprintf('LP1: %.3f [%.3f,%.3f]\nLP2: %.3f [%.3f,%.3f]', ...
            roc_data.LP1.AUC_mean, roc_data.LP1.AUC_CI(1), roc_data.LP1.AUC_CI(2), ...
            roc_data.LP2.AUC_mean, roc_data.LP2.AUC_CI(1), roc_data.LP2.AUC_CI(2));
        
        text(0.6, 0.2, auc_text, 'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'black');
        
        % Formatting
        xlabel('False Positive Rate','FontSize',10);
        ylabel('True Positive Rate','FontSize',10);
        title(sprintf('EDSS-%s', metric),'FontSize',12,'FontWeight','bold');
        axis square; xlim([0 1]); ylim([0 1]);
        set(gca,'FontSize',9,'LineWidth',1,'TickDir','out','Box','off');
        grid on; grid minor;
    end
end

% Row 2: LP11 and LP22 for EDSS
for metric_idx = 1:length(imaging_metrics)
    metric = imaging_metrics{metric_idx};
    field_name = sprintf('EDSS_%s', metric);
    
    if isfield(figure_data, field_name)
        roc_data = figure_data.(field_name);
        subplot(4, 5, 5 + metric_idx);
        hold on;
        
        % Extract data
        y = roc_data.outcome1;  % Use outcome1 for LP11
        prob_LP11 = roc_data.LP11.probabilities;
        prob_LP22 = roc_data.LP22.probabilities;
        
        % Plot ROC curves with confidence intervals
        % LP11 curve
        rm_lp11 = rocmetrics(y, prob_LP11, 1, 'NumBootstraps', 1000);
        h_lp11 = plot(rm_lp11, 'ShowConfidenceIntervals', true, 'Marker', 'none');
        h_lp11.LineWidth = 2;
        h_lp11.Color = color_lp11;
        h_lp11.LineStyle = '-';
        
        % LP22 curve
        rm_lp22 = rocmetrics(y, prob_LP22, 1, 'NumBootstraps', 1000);
        h_lp22 = plot(rm_lp22, 'ShowConfidenceIntervals', true, 'Marker', 'none');
        h_lp22.LineWidth = 2;
        h_lp22.Color = color_lp22;
        h_lp22.LineStyle = '-';
        
        % Chance diagonal
        plot([0 1],[0 1],'k--','LineWidth',1.5);
        
        % Add AUC (95% CI) text on the plot
        auc_text = sprintf('LP11: %.3f [%.3f,%.3f]\nLP22: %.3f [%.3f,%.3f]', ...
            roc_data.LP11.AUC_mean, roc_data.LP11.AUC_CI(1), roc_data.LP11.AUC_CI(2), ...
            roc_data.LP22.AUC_mean, roc_data.LP22.AUC_CI(1), roc_data.LP22.AUC_CI(2));
        
        text(0.6, 0.2, auc_text, 'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'black');
        
        % Formatting
        xlabel('False Positive Rate','FontSize',10);
        ylabel('True Positive Rate','FontSize',10);
        title(sprintf('EDSS-%s (Demo)', metric),'FontSize',12,'FontWeight','bold');
        axis square; xlim([0 1]); ylim([0 1]);
        set(gca,'FontSize',9,'LineWidth',1,'TickDir','out','Box','off');
        grid on; grid minor;
    end
end

% Row 3: LP1 and LP2 for MSPro
for metric_idx = 1:length(imaging_metrics)
    metric = imaging_metrics{metric_idx};
    field_name = sprintf('MSPro_%s', metric);
    
    if isfield(figure_data, field_name)
        roc_data = figure_data.(field_name);
        subplot(4, 5, 10 + metric_idx);
        hold on;
        
        % Extract data
        y = roc_data.outcome1;  % Use outcome1 for LP1
        prob_LP1 = roc_data.LP1.probabilities;
        prob_LP2 = roc_data.LP2.probabilities;
        
        % Plot ROC curves with confidence intervals
        % LP1 curve
        rm_lp1 = rocmetrics(y, prob_LP1, 1, 'NumBootstraps', 1000);
        h_lp1 = plot(rm_lp1, 'ShowConfidenceIntervals', true, 'Marker', 'none');
        h_lp1.LineWidth = 2;
        h_lp1.Color = color_lp1;
        h_lp1.LineStyle = '-';
        
        % LP2 curve
        rm_lp2 = rocmetrics(y, prob_LP2, 1, 'NumBootstraps', 1000);
        h_lp2 = plot(rm_lp2, 'ShowConfidenceIntervals', true, 'Marker', 'none');
        h_lp2.LineWidth = 2;
        h_lp2.Color = color_lp2;
        h_lp2.LineStyle = '-';
        
        % Chance diagonal
        plot([0 1],[0 1],'k--','LineWidth',1.5);
        
        % Add AUC (95% CI) text on the plot
        auc_text = sprintf('LP1: %.3f [%.3f,%.3f]\nLP2: %.3f [%.3f,%.3f]', ...
            roc_data.LP1.AUC_mean, roc_data.LP1.AUC_CI(1), roc_data.LP1.AUC_CI(2), ...
            roc_data.LP2.AUC_mean, roc_data.LP2.AUC_CI(1), roc_data.LP2.AUC_CI(2));
        
        text(0.6, 0.2, auc_text, 'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'black');
        
        % Formatting
        xlabel('False Positive Rate','FontSize',10);
        ylabel('True Positive Rate','FontSize',10);
        title(sprintf('MSPro-%s', metric),'FontSize',12,'FontWeight','bold');
        axis square; xlim([0 1]); ylim([0 1]);
        set(gca,'FontSize',9,'LineWidth',1,'TickDir','out','Box','off');
        grid on; grid minor;
    end
end

% Row 4: LP11 and LP22 for MSPro
for metric_idx = 1:length(imaging_metrics)
    metric = imaging_metrics{metric_idx};
    field_name = sprintf('MSPro_%s', metric);
    
    if isfield(figure_data, field_name)
        roc_data = figure_data.(field_name);
        subplot(4, 5, 15 + metric_idx);
        hold on;
        
        % Extract data
        y = roc_data.outcome1;  % Use outcome1 for LP11
        prob_LP11 = roc_data.LP11.probabilities;
        prob_LP22 = roc_data.LP22.probabilities;
        
        % Plot ROC curves with confidence intervals
        % LP11 curve
        rm_lp11 = rocmetrics(y, prob_LP11, 1, 'NumBootstraps', 1000);
        h_lp11 = plot(rm_lp11, 'ShowConfidenceIntervals', true, 'Marker', 'none');
        h_lp11.LineWidth = 2;
        h_lp11.Color = color_lp11;
        h_lp11.LineStyle = '-';
        
        % LP22 curve
        rm_lp22 = rocmetrics(y, prob_LP22, 1, 'NumBootstraps', 1000);
        h_lp22 = plot(rm_lp22, 'ShowConfidenceIntervals', true, 'Marker', 'none');
        h_lp22.LineWidth = 2;
        h_lp22.Color = color_lp22;
        h_lp22.LineStyle = '-';
        
        % Chance diagonal
        plot([0 1],[0 1],'k--','LineWidth',1.5);
        
        % Add AUC (95% CI) text on the plot
        auc_text = sprintf('LP11: %.3f [%.3f,%.3f]\nLP22: %.3f [%.3f,%.3f]', ...
            roc_data.LP11.AUC_mean, roc_data.LP11.AUC_CI(1), roc_data.LP11.AUC_CI(2), ...
            roc_data.LP22.AUC_mean, roc_data.LP22.AUC_CI(1), roc_data.LP22.AUC_CI(2));
        
        text(0.6, 0.2, auc_text, 'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'black');
        
        % Formatting
        xlabel('False Positive Rate','FontSize',10);
        ylabel('True Positive Rate','FontSize',10);
        title(sprintf('MSPro-%s (Demo)', metric),'FontSize',12,'FontWeight','bold');
        axis square; xlim([0 1]); ylim([0 1]);
        set(gca,'FontSize',9,'LineWidth',1,'TickDir','out','Box','off');
        grid on; grid minor;
    end
end

% Overall title
sgtitle('ROC Analysis - All Clinical Measures and Imaging Metrics (DEV VERSION)', 'FontSize', 16, 'FontWeight', 'bold');


end
