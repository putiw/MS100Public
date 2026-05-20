function paper_bh_p_correction()
% PAPER_BH_P_CORRECTION Performs Bonferroni-Holm correction across multiple Excel files
%
% This script:
% 1. Loads p-values from "Alternative_Layout" sheet in two Excel files
% 2. Performs Bonferroni-Holm correction across all p-values
% 3. Saves results to new Excel file with separate sheets for each model type
%
% Input files:
%   - paper_multi_bi_ridge.xlsx (Sheet: Alternative_Layout) - Binary/logistic models
%   - paper_multi_continuous_ridge.xlsx (Sheet: Alternative_Layout) - Continuous/linear models
%
% Output file:
%   - paper_multi_corrected.xlsx with sheets:
%     - Binary_Models
%     - Continuous_Models

clearvars; clc;

% Define input files (relative to repo root)
repo_root = fileparts(mfilename('fullpath'));
file1 = fullfile(repo_root, 'outputs', 'paper_multi_bi_ridge.xlsx');
file2 = fullfile(repo_root, 'outputs', 'paper_multi_continuous_ridge.xlsx');
sheet_name = 'Alternative_Layout';

fprintf('\n=== BONFERRONI-HOLM CORRECTION ===\n');
fprintf('Processing files:\n');
fprintf('  1. %s\n', file1);
fprintf('  2. %s\n', file2);
fprintf('  Sheet: %s\n\n', sheet_name);

% Check if files exist
if ~exist(file1, 'file')
    error('File not found: %s', file1);
end
if ~exist(file2, 'file')
    error('File not found: %s', file2);
end

% Read data from both files
fprintf('Reading binary/logistic models...\n');
data_binary = readtable(file1, 'Sheet', sheet_name);
fprintf('  Size: %d rows x %d columns\n', height(data_binary), width(data_binary));

fprintf('Reading continuous/linear models...\n');
data_continuous = readtable(file2, 'Sheet', sheet_name);
fprintf('  Size: %d rows x %d columns\n', height(data_continuous), width(data_continuous));

% Find p-value columns in both datasets
fprintf('\n=== IDENTIFYING P-VALUE COLUMNS ===\n');

var_names_binary = data_binary.Properties.VariableNames;
var_names_continuous = data_continuous.Properties.VariableNames;

p_value_cols_binary = {};
p_value_cols_continuous = {};

for i = 1:length(var_names_binary)
    if contains(lower(var_names_binary{i}), 'p_value') || contains(lower(var_names_binary{i}), 'p-value')
        p_value_cols_binary{end+1} = var_names_binary{i};
    end
end

for i = 1:length(var_names_continuous)
    if contains(lower(var_names_continuous{i}), 'p_value') || contains(lower(var_names_continuous{i}), 'p-value')
        p_value_cols_continuous{end+1} = var_names_continuous{i};
    end
end

fprintf('Binary models - p-value columns: %s\n', strjoin(p_value_cols_binary, ', '));
fprintf('Continuous models - p-value columns: %s\n', strjoin(p_value_cols_continuous, ', '));

% Collect all p-values from both files
fprintf('\n=== COLLECTING ALL P-VALUES ===\n');

all_p_values = [];
p_value_map = []; % [file_idx, col_idx, row_idx]

% Collect from binary models
for col_idx = 1:length(p_value_cols_binary)
    col_name = p_value_cols_binary{col_idx};
    p_vals = data_binary.(col_name);

    % Convert to numeric if needed
    if iscell(p_vals)
        p_vals_numeric = nan(size(p_vals));
        for i = 1:length(p_vals)
            if ischar(p_vals{i}) || isstring(p_vals{i})
                p_vals_numeric(i) = str2double(p_vals{i});
            elseif isnumeric(p_vals{i})
                p_vals_numeric(i) = p_vals{i};
            end
        end
        p_vals = p_vals_numeric;
    end

    % Find valid p-values (numeric, not NaN, between 0 and 1)
    valid_idx = ~isnan(p_vals) & p_vals >= 0 & p_vals <= 1;

    for row_idx = find(valid_idx)'
        all_p_values(end+1) = p_vals(row_idx);
        p_value_map(end+1, :) = [1, col_idx, row_idx]; % file 1 = binary
    end
end

n_binary_pvals = length(all_p_values);
fprintf('Binary models: %d valid p-values\n', n_binary_pvals);

% Collect from continuous models
for col_idx = 1:length(p_value_cols_continuous)
    col_name = p_value_cols_continuous{col_idx};
    p_vals = data_continuous.(col_name);

    % Convert to numeric if needed
    if iscell(p_vals)
        p_vals_numeric = nan(size(p_vals));
        for i = 1:length(p_vals)
            if ischar(p_vals{i}) || isstring(p_vals{i})
                p_vals_numeric(i) = str2double(p_vals{i});
            elseif isnumeric(p_vals{i})
                p_vals_numeric(i) = p_vals{i};
            end
        end
        p_vals = p_vals_numeric;
    end

    % Find valid p-values (numeric, not NaN, between 0 and 1)
    valid_idx = ~isnan(p_vals) & p_vals >= 0 & p_vals <= 1;

    for row_idx = find(valid_idx)'
        all_p_values(end+1) = p_vals(row_idx);
        p_value_map(end+1, :) = [2, col_idx, row_idx]; % file 2 = continuous
    end
end

n_continuous_pvals = length(all_p_values) - n_binary_pvals;
fprintf('Continuous models: %d valid p-values\n', n_continuous_pvals);
fprintf('TOTAL p-values for correction: %d\n', length(all_p_values));

% Perform corrections
fprintf('\n=== PERFORMING CORRECTIONS ===\n');

% Bonferroni-Holm correction
fprintf('Bonferroni-Holm correction...\n');
[p_values_holm, ~] = bonf_holm(all_p_values, 0.05);

% Report statistics
fprintf('\n=== CORRECTION STATISTICS ===\n');
alpha = 0.05;
n_sig_uncorrected = sum(all_p_values < alpha);
n_sig_holm = sum(p_values_holm < alpha);

fprintf('Results at α = %.2f:\n', alpha);
fprintf('  Uncorrected significant:      %2d / %d (%.1f%%)\n', ...
    n_sig_uncorrected, length(all_p_values), 100*n_sig_uncorrected/length(all_p_values));
fprintf('  Bonferroni-Holm significant:  %2d / %d (%.1f%%)\n', ...
    n_sig_holm, length(all_p_values), 100*n_sig_holm/length(all_p_values));

% Add corrected values back to tables
fprintf('\n=== ADDING CORRECTED VALUES TO TABLES ===\n');

% Process binary models
for col_idx = 1:length(p_value_cols_binary)
    col_name = p_value_cols_binary{col_idx};

    % Create new column name for Bonferroni-Holm correction
    % Use proper formatting for publication table
    holm_col_name = strrep(col_name, 'p_value', 'p_adj');
    holm_col_name = strrep(holm_col_name, 'p-value', 'p-adj');

    % Initialize with NaN
    holm_vals = nan(height(data_binary), 1);

    % Fill in corrected values
    for i = 1:size(p_value_map, 1)
        if p_value_map(i, 1) == 1 && p_value_map(i, 2) == col_idx
            row_idx = p_value_map(i, 3);
            holm_vals(row_idx) = p_values_holm(i);
        end
    end

    % Add column after original p-value column
    data_binary = addvars(data_binary, holm_vals, 'After', col_name, 'NewVariableNames', holm_col_name);

    fprintf('Binary models: Added %s\n', holm_col_name);
end

% Process continuous models
for col_idx = 1:length(p_value_cols_continuous)
    col_name = p_value_cols_continuous{col_idx};

    % Create new column name for Bonferroni-Holm correction
    % Use proper formatting for publication table
    holm_col_name = strrep(col_name, 'p_value', 'p_adj');
    holm_col_name = strrep(holm_col_name, 'p-value', 'p-adj');

    % Initialize with NaN
    holm_vals = nan(height(data_continuous), 1);

    % Fill in corrected values
    for i = 1:size(p_value_map, 1)
        if p_value_map(i, 1) == 2 && p_value_map(i, 2) == col_idx
            row_idx = p_value_map(i, 3);
            holm_vals(row_idx) = p_values_holm(i);
        end
    end

    % Add column after original p-value column
    data_continuous = addvars(data_continuous, holm_vals, 'After', col_name, 'NewVariableNames', holm_col_name);

    fprintf('Continuous models: Added %s\n', holm_col_name);
end

% Save to new Excel file
output_file = fullfile(repo_root, 'outputs', 'paper_multi_corrected.xlsx');

fprintf('\n=== SAVING RESULTS ===\n');
fprintf('Output file: %s\n', output_file);

writetable(data_binary, output_file, 'Sheet', 'Binary_Models');
fprintf('Saved sheet: Binary_Models\n');

writetable(data_continuous, output_file, 'Sheet', 'Continuous_Models');
fprintf('Saved sheet: Continuous_Models\n');

% Create summary sheet for EDSS and MSPro
fprintf('\n=== CREATING SUMMARY SHEET ===\n');
summary_table = createSummaryTable(data_binary);
writetable(summary_table, output_file, 'Sheet', 'EDSS_MSPro_Summary');
fprintf('Saved sheet: EDSS_MSPro_Summary\n');

fprintf('\n=== CORRECTION COMPLETE ===\n');
fprintf('Corrected p-values saved to: %s\n', output_file);
fprintf('Total tests corrected: %d (%d binary + %d continuous)\n', ...
    length(all_p_values), n_binary_pvals, n_continuous_pvals);
fprintf('\nNew column added:\n');
fprintf('  - p_adj: Bonferroni-Holm adjusted p-values\n');
fprintf('\nNote: Bonferroni-Holm correction controls family-wise error rate\n');
fprintf('      while being more powerful than standard Bonferroni correction.\n');
fprintf('\nFor table footnote, use:\n');
fprintf('  "p_adj: p-values adjusted for multiple comparisons using the\n');
fprintf('   Bonferroni-Holm method (Holm, 1979)"\n');

end

function summary_table = createSummaryTable(data_binary)
    % Create summary table with EDSS and MSPro results
    % Rows: EDSS-Lnorm, EDSS-T1, EDSS-MTR, EDSS-FA, EDSS-MD, MSPro-Lnorm, MSPro-T1, MSPro-MTR, MSPro-FA, MSPro-MD
    % Columns: Tract-based AUC (95% CI), Classic-based AUC (95% CI), Tract-based AIC, Classic-based AIC

    % Define the combinations
    measures = {'EDSS', 'EDSS', 'EDSS', 'EDSS', 'EDSS', 'MSPro', 'MSPro', 'MSPro', 'MSPro', 'MSPro'};
    metrics = {'Lnorm', 'T1', 'MTR', 'FA', 'MD', 'Lnorm', 'T1', 'MTR', 'FA', 'MD'};

    % Initialize output arrays
    row_labels = cell(10, 1);
    tract_n = cell(10, 1);
    classic_n = cell(10, 1);
    tract_auc_ci = cell(10, 1);
    classic_auc_ci = cell(10, 1);
    tract_aic = nan(10, 1);
    classic_aic = nan(10, 1);

    % Find the actual column names for N, AIC and AUC_CI
    col_names = data_binary.Properties.VariableNames;
    n_col = '';
    aic_col = '';
    auc_ci_col = '';

    for j = 1:length(col_names)
        col_lower = lower(col_names{j});
        if strcmp(col_lower, 'n')
            n_col = col_names{j};
        end
        if contains(col_lower, 'aic')
            aic_col = col_names{j};
        end
        if contains(col_lower, 'auc')
            auc_ci_col = col_names{j};
        end
    end

    fprintf('Using columns: N=%s, AIC=%s, AUC_CI=%s\n', n_col, aic_col, auc_ci_col);

    if isempty(auc_ci_col)
        error('Could not find AUC_CI column in data');
    end

    % Extract data for each combination
    for i = 1:10
        measure = measures{i};
        metric = metrics{i};
        row_labels{i} = sprintf('%s-%s', measure, metric);

        % Find matching rows in data_binary
        clinical_metric = data_binary.ClinicalMetric;

        % Handle both string and cell array formats
        if iscell(clinical_metric)
            match_str = sprintf('%s-%s', measure, metric);
            tract_idx = find(strcmp(clinical_metric, match_str) & strcmp(data_binary.Model, 'tract-based'));
            classic_idx = find(strcmp(clinical_metric, match_str) & strcmp(data_binary.Model, 'classic'));
        else
            match_str = sprintf('%s-%s', measure, metric);
            tract_idx = find(contains(string(clinical_metric), match_str) & contains(string(data_binary.Model), 'tract-based'));
            classic_idx = find(contains(string(clinical_metric), match_str) & contains(string(data_binary.Model), 'classic'));
        end

        % Extract tract-based values
        if ~isempty(tract_idx)
            % N (sample size)
            if ~isempty(n_col)
                n_val = data_binary.(n_col)(tract_idx(1));
                if iscell(n_val)
                    tract_n{i} = n_val{1};
                else
                    tract_n{i} = num2str(n_val);
                end
            else
                tract_n{i} = 'N/A';
            end

            % AUC_CI column contains the full string like "0.750 [0.650, 0.850]"
            auc_ci_str = data_binary.(auc_ci_col)(tract_idx(1));

            % Convert cell to string if needed
            if iscell(auc_ci_str)
                auc_ci_str = auc_ci_str{1};
            end

            % Format it properly
            tract_auc_ci{i} = formatAucCI(auc_ci_str);

            % AIC
            if ~isempty(aic_col)
                tract_aic(i) = extractNumeric(data_binary.(aic_col)(tract_idx(1)));
            end
        else
            tract_n{i} = 'N/A';
            tract_auc_ci{i} = 'N/A';
        end

        % Extract classic-based values
        if ~isempty(classic_idx)
            % N (sample size)
            if ~isempty(n_col)
                n_val = data_binary.(n_col)(classic_idx(1));
                if iscell(n_val)
                    classic_n{i} = n_val{1};
                else
                    classic_n{i} = num2str(n_val);
                end
            else
                classic_n{i} = 'N/A';
            end

            % AUC_CI column contains the full string like "0.750 [0.650, 0.850]"
            auc_ci_str = data_binary.(auc_ci_col)(classic_idx(1));

            % Convert cell to string if needed
            if iscell(auc_ci_str)
                auc_ci_str = auc_ci_str{1};
            end

            % Format it properly
            classic_auc_ci{i} = formatAucCI(auc_ci_str);

            % AIC
            if ~isempty(aic_col)
                classic_aic(i) = extractNumeric(data_binary.(aic_col)(classic_idx(1)));
            end
        else
            classic_n{i} = 'N/A';
            classic_auc_ci{i} = 'N/A';
        end
    end

    % Create table
    summary_table = table(row_labels, tract_n, classic_n, tract_auc_ci, classic_auc_ci, tract_aic, classic_aic, ...
        'VariableNames', {'Clinical_Metric', 'TractBased_N', 'ClassicBased_N', 'TractBased_AUC_95CI', 'ClassicBased_AUC_95CI', 'TractBased_AIC', 'ClassicBased_AIC'});

    fprintf('Summary table created with %d rows\n', height(summary_table));
end

function formatted_str = formatAucCI(auc_ci_str)
    % Format AUC with CI string from formats like:
    % "0.750 [0.650, 0.850]" -> "0.750 (0.650-0.850)"
    % "0.750 [0.650,0.850]" -> "0.750 (0.650-0.850)"

    if isempty(auc_ci_str) || strcmp(auc_ci_str, 'N/A')
        formatted_str = 'N/A';
        return;
    end

    % Convert to string
    str_val = char(auc_ci_str);

    % Try to extract AUC value and CI
    % Pattern: "AUC [CI_low, CI_high]"
    tokens = regexp(str_val, '([\d.]+)\s*\[([\d.]+)\s*,\s*([\d.]+)\]', 'tokens');

    if ~isempty(tokens)
        auc = tokens{1}{1};
        ci_low = tokens{1}{2};
        ci_high = tokens{1}{3};
        formatted_str = sprintf('%s (%s-%s)', auc, ci_low, ci_high);
    else
        % If parsing fails, return as is
        formatted_str = str_val;
    end
end

function val = extractNumeric(cell_val)
    % Extract numeric value from cell or string
    if isnumeric(cell_val)
        val = cell_val;
    elseif iscell(cell_val)
        val = str2double(cell_val{1});
    elseif ischar(cell_val) || isstring(cell_val)
        val = str2double(cell_val);
    else
        val = NaN;
    end
end

