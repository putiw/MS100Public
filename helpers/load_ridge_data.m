function [merged_data, preds1, preds2, data1, data2, y_raw, Demo, clinical_measure, demographic_vars] = load_ridge_data(varargin)
% LOAD_RIDGE_DATA Load and prepare data for ridge regression analysis
%
% Usage:
%   [merged_data, preds1, preds2, data1, data2, y_raw, Demo, clinical_measure, demographic_vars] = load_ridge_data()
%   [merged_data, preds1, preds2, data1, data2, y_raw, Demo, clinical_measure, demographic_vars] = load_ridge_data('measure', 'EDSS', 'model1_metric', 'MTR')
%
% Parameters:
%   'measure'        - Clinical measure ('EDSS', 'T25FW', 'x9HPTD', 'x9HPTND', 'MSPro')
%   'model1_metric' - NRF metric ('MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm')
%   'model2_metric'  - CRF metric ('MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm')
%   'demographics'   - Demographic variables ({'Age', 'Gender', 'DurationOfDisease'})
%   'check_predictors' - Which predictors to check for missing data ('both', 'model1', 'model2', 'none')
%
% Returns:
%   merged_data      - Cleaned and merged data table
%   preds1           - Cell array of Model 1 predictor variable names
%   preds2           - Cell array of Model 2 predictor variable names
%   data1            - Matrix of Model 1 predictor data
%   data2            - Matrix of Model 2 predictor data
%   y_raw            - Raw clinical measure values
%   Demo             - Matrix of demographic covariates
%   clinical_measure - Clinical measure name
%   demographic_vars - Demographic variable names

% Parse input arguments
p = inputParser;
addParameter(p, 'measure', 'T25FW', @(x) ismember(x, {'EDSS', 'T25FW', 'x9HPTD', 'x9HPTND', 'MSPro'}));
addParameter(p, 'model1_metric', 'T1', @(x) ismember(x, {'MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm'}));
addParameter(p, 'model2_metric', 'T1', @(x) ismember(x, {'MTR', 'T1', 'FA', 'MD', 'AD', 'RD', 'LN', 'LV', 'Lnorm'}));
addParameter(p, 'demographics', {'Age', 'Gender'}, @iscell); %, 'DurationOfDisease'
addParameter(p, 'check_predictors', 'both', @(x) ismember(x, {'both', 'model1', 'model2', 'none'}));
parse(p, varargin{:});

% Extract parameters
clinical_measure = p.Results.measure;
model1_metric = p.Results.model1_metric;
model2_metric = p.Results.model2_metric;
demographic_vars = p.Results.demographics;
check_predictors = p.Results.check_predictors;

% Add helper functions
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'helpers')));

% Load configuration
config = load_config();
stats_dir = fullfile(config.bidsDir, config.statsDir);

% Load clinical data
clinical_data = readtable(fullfile(stats_dir, 'clinicalScore.xlsx'));

% Check if clinical measure exists in data
if ~ismember(clinical_measure, clinical_data.Properties.VariableNames)
    error('Clinical measure "%s" not found in clinical data. Available measures: %s', ...
        clinical_measure, strjoin(clinical_data.Properties.VariableNames, ', '));
end

% Define tract groups for Model 1
tract_groups = {
    'Association', ...
    'Cerebellar', ...
    'Occipitoparietal', ...
    'PB'  % ProjectionBrain (shortened)
};

% Define Icometrix metrics for Model 2
ico_metrics = {
    'periventricular', ...
    'juxtacortical', ...
    'infratentorial', ...
    'deepwhitematter'
};

% Load Model 1 data based on metric type
if ismember(model1_metric, {'LN', 'LV', 'Lnorm'})
    % Load lesion load data
    tract_data = readtable(fullfile(stats_dir, 'GroupTractLesionLoad.xlsx'));
else
    % Load stats data
    tract_data = readtable(fullfile(stats_dir, ['GroupTract' model1_metric '_All.xlsx']));
end

% Load Model 2 data based on metric type
if ismember(model2_metric, {'LN', 'LV', 'Lnorm'})
    % Load lesion load data
    ico_data = readtable(fullfile(stats_dir, 'icometrixLesionLoad.xlsx'));
else
    % Load stats data - try different possible filename patterns
    possible_files = {
        ['icometrix' model2_metric '.xlsx'],
        ['icometrixLesion' model2_metric '.xlsx'],
        'icometrixLesionStats.xlsx'
    };

    ico_data = [];
    for i = 1:length(possible_files)
        filepath = fullfile(stats_dir, possible_files{i});
        if exist(filepath, 'file')
            ico_data = readtable(filepath);
            fprintf('Loaded CRF data from: %s\n', possible_files{i});
            break;
        end
    end

    if isempty(ico_data)
        error('Could not find icometrix data file for metric: %s. Tried: %s', ...
            model2_metric, strjoin(possible_files, ', '));
    end
end

% Merge tables on SubjectID
merged_data = outerjoin(clinical_data, tract_data, 'Keys','SubjectID','MergeKeys',true);
merged_data = outerjoin(merged_data, ico_data, 'Keys','SubjectID','MergeKeys',true);
merged_data.Properties.VariableNames{ ...
    strcmp(merged_data.Properties.VariableNames,'SubjectID')} = 'Subject';

% Exclude controls and unwanted subjects
if ismember(clinical_measure,{'EDSS','MSPro'})
    merged_data = merged_data(~startsWith(merged_data.Subject,'sub-C'), :);
end 


% Construct predictor variables for Model 1
if ismember(model1_metric, {'LN', 'LV', 'Lnorm'})
    % For Lesion Load, construct variable names based on metric
    preds1 = {};
    for i = 1:length(tract_groups)
        group = tract_groups{i};
        preds1{end+1} = [group 'L' model1_metric];  % Left
        preds1{end+1} = [group 'R' model1_metric];  % Right
    end
else
    % For stats metrics, construct variable names based on tissue
    tract_groups_full = {
        'Association', ...
        'Cerebellar', ...
        'Occipitoparietal', ...
        'ProjectionBrainstem'  % Full name for stats
    };
    preds1 = {};
    for i = 1:length(tract_groups_full)
        group = tract_groups_full{i};
        preds1{end+1} = [group 'L_Tail'];
        preds1{end+1} = [group 'R_Tail'];
        preds1{end+1} = [group 'L_NAWM'];
        preds1{end+1} = [group 'R_NAWM'];
    end
end

% Construct predictor variables for Model 2
preds2 = {};
for i = 1:length(ico_metrics)
    preds2{end+1} = [ico_metrics{i} model2_metric];
end

% Validate that all predictor variables exist
missing_vars1 = setdiff(preds1, merged_data.Properties.VariableNames);
if ~isempty(missing_vars1)
    warning('Missing Model 1 predictor variables: %s', strjoin(missing_vars1, ', '));
    preds1 = setdiff(preds1, missing_vars1);
end

missing_vars2 = setdiff(preds2, merged_data.Properties.VariableNames);
if ~isempty(missing_vars2)
    warning('Missing Model 2 predictor variables: %s', strjoin(missing_vars2, ', '));
    preds2 = setdiff(preds2, missing_vars2);
end

% Check for NaN and zero values in selected variables
% Build list of variables to check based on check_predictors parameter
vars_to_check = [{clinical_measure}, demographic_vars];

% Track which predictors should have zero-checking
% Intensity metrics (MTR, T1, FA, MD, AD, RD) should exclude zeros
% Lesion metrics (LN, LV, Lnorm) should allow zeros (no lesions is valid)
intensity_metrics = {'MTR', 'T1', 'FA', 'MD', 'AD', 'RD'};
preds_to_check_zero1 = {};
preds_to_check_zero2 = {};

switch check_predictors
    case 'both'
        vars_to_check = [vars_to_check, preds1, preds2];
        if ismember(model1_metric, intensity_metrics)
            preds_to_check_zero1 = preds1;
        end
        if ismember(model2_metric, intensity_metrics)
            preds_to_check_zero2 = preds2;
        end
    case 'model1'
        vars_to_check = [vars_to_check, preds1];
        if ismember(model1_metric, intensity_metrics)
            preds_to_check_zero1 = preds1;
        end
    case 'model2'
        vars_to_check = [vars_to_check, preds2];
        if ismember(model2_metric, intensity_metrics)
            preds_to_check_zero2 = preds2;
        end
    case 'none'
        % Only check clinical measure and demographics
end

% Combine all predictors that need zero-checking
preds_to_check_zero = [preds_to_check_zero1, preds_to_check_zero2];

% Create exclusion mask for NaN and zero values
mask = false(height(merged_data),1);
for ii = 1:numel(vars_to_check)
    v = vars_to_check{ii};
    if ismember(v, merged_data.Properties.VariableNames)
        var_data = merged_data.(v);
        if isnumeric(var_data)
            % Always exclude NaN
            mask = mask | isnan(var_data);

            % For intensity metrics, also exclude zeros
            if ismember(v, preds_to_check_zero)
                mask = mask | (var_data == 0);
            end
        elseif iscell(var_data)
            mask = mask | cellfun(@(x) isempty(x) || (isnumeric(x) && isnan(x)), var_data);
        elseif iscategorical(var_data)
            mask = mask | isundefined(var_data);
        end
    end
end

% Report exclusions
n_excluded = sum(mask);
if n_excluded > 0
    fprintf('Excluding %d subjects due to missing/zero data in checked variables\n', n_excluded);
end

merged_data(mask, :) = [];

% Get raw clinical measure values
y_raw = merged_data.(clinical_measure);

% Get data for each model
data1 = merged_data{:, preds1};
data2 = merged_data{:, preds2};

% Prepare demographic covariates
Demo = [];
if ~isempty(demographic_vars)
    for i = 1:length(demographic_vars)
        var = demographic_vars{i};
        if ~ismember(var, merged_data.Properties.VariableNames)
            warning('Demographic variable "%s" not found in data. Skipping.', var);
            continue;
        end
        if strcmp(var, 'Gender')
            Demo = [Demo, double(strcmp(merged_data.Gender, 'M'))];  % 1 for Male, 0 for Female
        else
            Demo = [Demo, merged_data.(var)];
        end
    end
end

% Display data loading summary
fprintf('\n=== DATA LOADING SUMMARY ===\n');
fprintf('Clinical measure: %s\n', clinical_measure);
fprintf('Model 1 metric: %s (%d predictors)\n', model1_metric, length(preds1));
fprintf('Model 2 metric: %s (%d predictors)\n', model2_metric, length(preds2));
fprintf('Sample size after exclusions: %d\n', height(merged_data));
fprintf('Demographics: %s\n', strjoin(demographic_vars, ', '));
fprintf('=============================\n\n');

end
