function outputs = generate_manuscript_tables(varargin)
% GENERATE_MANUSCRIPT_TABLES Write the five manuscript tables as plain CSV.
%
% The table grids match Main Tables 1-2 and Supplementary Tables 3-5.
% No Word formatting, captions, or participant-level data are written.
%
% Minimal call:
%   generate_manuscript_tables( ...
%       'StatsDir', '/private/MsBIDS/derivatives/derivativesStats');

% Main Table 1 is calculated from clinicalScore.xlsx. Main Table 2 and
% Supplementary Tables 3-4 are formatted from the unified qMRI model CSV.
% Supplementary Table 5 is formatted from the aggregate lesion-load model CSV.


p = inputParser;
script_dir = fileparts(mfilename('fullpath'));
package_root = fileparts(script_dir);

addParameter(p, 'StatsDir', '', @is_text_scalar);
addParameter(p, 'ClinicalFile', '', @is_text_scalar);
addParameter(p, 'ModelFile', fullfile(package_root, 'results', ...
    'model_results', 'all_qmri_manuscript_models.csv'), @is_text_scalar);
addParameter(p, 'LesionModelFile', fullfile(package_root, 'results', ...
    'model_results', 'lesionload_manuscript_models.csv'), @is_text_scalar);
addParameter(p, 'OutputDir', fullfile(package_root, 'results', 'tables'), ...
    @is_text_scalar);
parse(p, varargin{:});

stats_dir = char(p.Results.StatsDir);
clinical_file = char(p.Results.ClinicalFile);
if isempty(clinical_file)
    assert(~isempty(stats_dir), ...
        'StatsDir or ClinicalFile is required to generate Main Table 1.');
    clinical_file = fullfile(stats_dir, 'clinicalScore.xlsx');
end
model_file = require_file(char(p.Results.ModelFile), 'all-qMRI model table');
lesion_file = require_file(char(p.Results.LesionModelFile), ...
    'lesion-load model table');
clinical_file = require_file(clinical_file, 'clinical score workbook');
output_dir = char(p.Results.OutputDir);
if ~isfolder(output_dir)
    mkdir(output_dir);
end

clinical = readtable(clinical_file, 'VariableNamingRule', 'preserve');
clinical = add_msfc_sdmt(clinical);
models = readtable(model_file, 'VariableNamingRule', 'preserve', ...
    'TextType', 'string');
lesion_models = readtable(lesion_file, 'VariableNamingRule', 'preserve', ...
    'TextType', 'string');

validate_qmri_models(models);
validate_lesion_models(lesion_models);

tables = struct();
tables.Main_Table_1 = build_main_table_1(clinical);
tables.Main_Table_2 = build_binary_performance_table(models, false);
tables.Supplementary_Table_3 = build_binary_performance_table(models, true);
tables.Supplementary_Table_4 = build_adjusted_p_table(models);
tables.Supplementary_Table_5 = build_lesion_performance_table(lesion_models);

expected = struct( ...
    'Main_Table_1', [13 5], ...
    'Main_Table_2', [15 5], ...
    'Supplementary_Table_3', [15 5], ...
    'Supplementary_Table_4', [99 6], ...
    'Supplementary_Table_5', [58 6]);

names = fieldnames(tables);
outputs = struct();
for i = 1:numel(names)
    name = names{i};
    values = tables.(name);
    assert(isequal(size(values), expected.(name)), ...
        '%s has an unexpected grid size.', name);
    output_path = fullfile(output_dir, [name '.csv']);
    write_plain_csv(output_path, values);
    outputs.(name) = output_path;
    fprintf('%s\n', output_path);
end
end


function values = build_main_table_1(T)
required = {'SubjectID','Group','Age','Gender','DurationOfDisease','EDSS', ...
    'T25FW','x9HPTD','x9HPTND','SDMTcorrect','MSFC_SDMT','MSPro','DrugLine'};
assert(isempty(setdiff(required, T.Properties.VariableNames)), ...
    'clinicalScore.xlsx is missing required Table 1 fields.');

group = numeric_vector(T.Group);
group_codes = 0:3;
expected_n = [43 49 17 23];
group_names = {'Healthy Controls','RRMS','PPMS','SPMS'};
for i = 1:numel(group_codes)
    assert(sum(group == group_codes(i)) == expected_n(i), ...
        'Unexpected participant count for %s.', group_names{i});
end
assert(height(T) == sum(expected_n), ...
    'Main Table 1 requires the final 132-participant clinical workbook.');

values = cell(13,5);
values(1,:) = [{'Parameter'}, group_names];
row_names = {'N (% Female)','Age','Disease duration (years)','EDSS', ...
    'Therapy: 1st line','Therapy: 2nd line','T25FW (s)','9HPT-D (s)', ...
    '9HPT-ND (s)','SDMT correct responses','MSFC-SDMT (z)','MSPro'};
values(2:end,1) = row_names(:);

gender = upper(strtrim(string(T.Gender)));
for i = 1:numel(group_codes)
    keep = group == group_codes(i);
    n = sum(keep);
    female_n = sum(gender(keep) == "F");
    values{2,i+1} = sprintf('%d (%d%%)', n, round(100*female_n/n));
    values{3,i+1} = mean_sd(T.Age(keep));
    values{8,i+1} = mean_sd(T.T25FW(keep));
    values{9,i+1} = mean_sd(T.x9HPTD(keep));
    values{10,i+1} = mean_sd(T.x9HPTND(keep));
    values{11,i+1} = mean_sd(T.SDMTcorrect(keep));
    msfc = numeric_vector(T.MSFC_SDMT(keep));
    values{12,i+1} = sprintf('%s (n=%d)', mean_sd(msfc), sum(isfinite(msfc)));

    if group_codes(i) == 0
        values(4:7,i+1) = repmat({'-'}, 4, 1);
        values{13,i+1} = '-';
    else
        values{4,i+1} = mean_sd(T.DurationOfDisease(keep));
        values{5,i+1} = median_iqr(T.EDSS(keep));
        values{6,i+1} = therapy_summary(T.DrugLine(keep), 1, n);
        values{7,i+1} = therapy_summary(T.DrugLine(keep), 2, n);
        values{13,i+1} = median_iqr(T.MSPro(keep));
    end
end

msfc_counts = zeros(1,4);
for i = 1:4
    msfc_counts(i) = sum(isfinite(numeric_vector(T.MSFC_SDMT(group == group_codes(i)))));
end
assert(isequal(msfc_counts, [43 49 16 19]), ...
    'Unexpected group-specific MSFC-SDMT complete-case counts.');
end


function values = build_binary_performance_table(T, adjusted)
outcomes = {'EDSS','MSPro'};
stored_metrics = {'T1','MTR','FA','MD','NDI','ODI','FWF'};
display_metrics = {'T1','MTR','FA','MD','NDI','ODI','ISOVF'};

values = cell(15,5);
if adjusted
    values(1,:) = {'Clinical Metric (n=number of subjects)', ...
        'Tract-Based AUC (95% CI)', 'Classical-Based AUC (95% CI)', ...
        'Tract-Based AIC', 'Classical-Based AIC'};
    separator = char(8211);
else
    values(1,:) = {'Clinical Metric (n=number of patients)', ...
        'Tract-based AUC (95% CI)', 'Classical-based AUC (95% CI)', ...
        'Tract-based AIC', 'Classical-based AIC'};
    separator = '-';
end

row = 2;
for oi = 1:numel(outcomes)
    for mi = 1:numel(stored_metrics)
        tract = one_qmri_row(T, outcomes{oi}, stored_metrics{mi}, ...
            'Tract-based', adjusted);
        classic = one_qmri_row(T, outcomes{oi}, stored_metrics{mi}, ...
            'Classical-region', adjusted);
        n_tract = scalar_number(tract.N);
        n_classic = scalar_number(classic.N);
        p_tract = scalar_number(tract.PositiveN);
        p_classic = scalar_number(classic.PositiveN);
        assert(n_tract == n_classic && p_tract == p_classic, ...
            'Binary tract/classical case sets do not match.');

        values{row,1} = sprintf('%s-%s (n=%d)', outcomes{oi}, ...
            display_metrics{mi}, n_tract);
        values{row,2} = auc_text(tract, separator);
        values{row,3} = auc_text(classic, separator);
        values{row,4} = sprintf('%.2f', scalar_number(tract.AIC));
        values{row,5} = sprintf('%.2f', scalar_number(classic.AIC));
        row = row + 1;
    end
end
end


function values = build_adjusted_p_table(T)
outcomes = {'EDSS','MSPro','T25FW','9HPT-D','9HPT-ND','SDMT','MSFC-SDMT'};
stored_metrics = {'T1','MTR','FA','MD','NDI','ODI','FWF'};
display_metrics = {'T1','MTR','FA','MD','NDI','ODI','ISOVF'};
models = {'Tract-based','Classical-region'};

values = cell(99,6);
values(1,:) = {'Outcome','MRI Metric','Model','N','p','p_adj'};
row = 2;
for oi = 1:numel(outcomes)
    for mi = 1:numel(stored_metrics)
        for fi = 1:numel(models)
            result = one_qmri_row(T, outcomes{oi}, stored_metrics{mi}, ...
                models{fi}, true);
            if mi == 1 && fi == 1
                values{row,1} = outcomes{oi};
            else
                values{row,1} = '';
            end
            if fi == 1
                values{row,2} = display_metrics{mi};
                values{row,3} = 'Tract-based';
            else
                values{row,2} = '';
                values{row,3} = 'Classic';
            end
            n = scalar_number(result.N);
            if ismember(outcomes{oi}, {'EDSS','MSPro'})
                values{row,4} = sprintf('%d (%d)', n, ...
                    scalar_number(result.PositiveN));
            else
                values{row,4} = sprintf('%d', n);
            end
            raw_p = scalar_number(result.p_raw);
            if ismember(stored_metrics{mi}, {'T1','MTR','FA','MD'})
                values{row,5} = legacy_p_text(raw_p);
            else
                values{row,5} = p_text(raw_p);
            end
            values{row,6} = p_text(scalar_number(result.p_adj));
            row = row + 1;
        end
    end
end
end


function values = build_lesion_performance_table(T)
binary_outcomes = {'EDSS','MSPro'};
continuous_outcomes = {'T25FW','9HPT-D','9HPT-ND','SDMT','MSFC-SDMT'};
metrics = {'LN','LV','LN','LN','LV','LV','Lnorm','Lnorm'};
models = {'Whole Brain','Whole Brain','Tract-based','Classical-region', ...
    'Tract-based','Classical-region','Tract-based','Classical-region'};

values = cell(58,6);
values(1,:) = {'Outcome','Metric','Model','N','AUC (95% CI)','AIC'};
row = 2;
for oi = 1:numel(binary_outcomes)
    for i = 1:numel(metrics)
        result = one_lesion_row(T, binary_outcomes{oi}, metrics{i}, models{i});
        if i == 1
            values{row,1} = binary_outcomes{oi};
        else
            values{row,1} = '';
        end
        values{row,2} = metrics{i};
        values{row,3} = lesion_model_label(models{i});
        values{row,4} = sprintf('%d (%d)', scalar_number(result.N), ...
            scalar_number(result.PositiveN));
        values{row,5} = sprintf('%.3f [%.3f, %.3f]', ...
            scalar_number(result.Performance), scalar_number(result.CI_lower), ...
            scalar_number(result.CI_upper));
        values{row,6} = sprintf('%.2f', scalar_number(result.AIC));
        row = row + 1;
    end
end

values(row,:) = {'Outcome','Metric','Model','N','R2','AIC'};
row = row + 1;
for oi = 1:numel(continuous_outcomes)
    for i = 1:numel(metrics)
        result = one_lesion_row(T, continuous_outcomes{oi}, metrics{i}, models{i});
        if i == 1
            values{row,1} = continuous_outcomes{oi};
        else
            values{row,1} = '';
        end
        values{row,2} = metrics{i};
        values{row,3} = lesion_model_label(models{i});
        values{row,4} = sprintf('%d', scalar_number(result.N));
        values{row,5} = sprintf('%.3f', scalar_number(result.Performance));
        values{row,6} = sprintf('%.2f', scalar_number(result.AIC));
        row = row + 1;
    end
end
assert(row == 59, 'Supplementary Table 5 row construction failed.');
end


function validate_qmri_models(T)
required = {'Outcome','Metric','Model','Adjusted','Demographics','N', ...
    'PositiveN','Performance','CI_lower','CI_upper','p_raw','p_adj', ...
    'MultiplicityFamilyN','AIC'};
assert(isempty(setdiff(required, T.Properties.VariableNames)), ...
    'all_qmri_manuscript_models.csv has unexpected columns.');
assert(height(T) == 196, 'Expected 196 all-qMRI model rows.');

adjusted = numeric_vector(T.Adjusted);
assert(sum(adjusted == 1) == 98 && sum(adjusted == 0) == 98, ...
    'Expected 98 adjusted and 98 unadjusted qMRI rows.');
family_n = numeric_vector(T.MultiplicityFamilyN);
assert(all(family_n(adjusted == 1) == 98), ...
    'Adjusted qMRI rows must use the 98-test family.');

raw_p = numeric_vector(T.p_raw(adjusted == 1));
adjusted_p = numeric_vector(T.p_adj(adjusted == 1));
assert(max(abs(adjusted_p - min(98*raw_p,1))) < 1e-10, ...
    'Adjusted p values do not match Bonferroni m=98.');
assert(sum(adjusted_p < 0.05) == 46, ...
    'Unexpected number of significant adjusted qMRI rows.');

outcome = string(T.Outcome);
demographics = string(T.Demographics);
binary = ismember(outcome, ["EDSS","MSPro"]);
assert(all(demographics(adjusted == 1 & binary) == ...
    "Age+GenderNum+DurationOfDisease"), ...
    'Adjusted binary models must use age, gender, and disease duration.');
assert(all(demographics(adjusted == 1 & ~binary) == "Age+GenderNum"), ...
    'Adjusted continuous models must use age and gender.');

outcomes = ["EDSS","MSPro","T25FW","9HPT-D","9HPT-ND","SDMT","MSFC-SDMT"];
metrics = ["T1","MTR","FA","MD","NDI","ODI","FWF"];
models = ["Tract-based","Classical-region"];
for oi = 1:numel(outcomes)
    for mi = 1:numel(metrics)
        for fi = 1:numel(models)
            key = outcome == outcomes(oi) & string(T.Metric) == metrics(mi) & ...
                string(T.Model) == models(fi);
            assert(sum(key) == 2 && isequal(sort(adjusted(key)), [0;1]), ...
                'Missing or duplicate qMRI model key.');
        end
    end
end
end


function validate_lesion_models(T)
required = {'Outcome','Metric','Model','Adjusted','Demographics','N', ...
    'PositiveN','PerformanceName','Performance','CI_lower','CI_upper','p_raw','AIC'};
assert(isempty(setdiff(required, T.Properties.VariableNames)), ...
    'lesionload_manuscript_models.csv has unexpected columns.');
assert(height(T) == 56, 'Expected 56 lesion-load model rows.');
assert(all(numeric_vector(T.Adjusted) == 0), ...
    'Supplementary Table 5 uses unadjusted lesion-load models.');
assert(all(string(T.Demographics) == "none"), ...
    'Supplementary Table 5 lesion-load rows must be unadjusted.');

binary = ismember(string(T.Outcome), ["EDSS","MSPro"]);
assert(sum(binary) == 16 && sum(~binary) == 40, ...
    'Expected 16 binary and 40 continuous lesion-load rows.');
assert(all(string(T.PerformanceName(binary)) == "AUC") && ...
       all(string(T.PerformanceName(~binary)) == "R2"), ...
    'Unexpected lesion-load performance labels.');
end


function row = one_qmri_row(T, outcome, metric, model, adjusted)
keep = string(T.Outcome) == string(outcome) & ...
    string(T.Metric) == string(metric) & ...
    string(T.Model) == string(model) & ...
    numeric_vector(T.Adjusted) == double(adjusted);
assert(sum(keep) == 1, ...
    'Expected one qMRI row for %s/%s/%s/adjusted=%d.', ...
    outcome, metric, model, adjusted);
row = T(keep,:);
end


function row = one_lesion_row(T, outcome, metric, model)
keep = string(T.Outcome) == string(outcome) & ...
    string(T.Metric) == string(metric) & ...
    string(T.Model) == string(model);
assert(sum(keep) == 1, ...
    'Expected one lesion row for %s/%s/%s.', outcome, metric, model);
row = T(keep,:);
end


function value = auc_text(row, separator)
value = sprintf('%.3f (%.3f%s%.3f)', scalar_number(row.Performance), ...
    scalar_number(row.CI_lower), separator, scalar_number(row.CI_upper));
end


function value = lesion_model_label(model)
if strcmp(model, 'Classical-region')
    value = 'Classic';
else
    value = model;
end
end


function value = p_text(p_value)
if p_value < 0.0001
    value = '< 0.0001';
else
    value = sprintf('%.4f', p_value);
end
end


function value = legacy_p_text(p_value)
% The retained manuscript rows first rounded to four decimal places.
value = sprintf('%.4f', p_value);
if strcmp(value, '0.0000')
    value = '< 0.0001';
end
end


function value = mean_sd(values)
values = numeric_vector(values);
values = values(isfinite(values));
assert(~isempty(values), 'Cannot summarize an empty clinical field.');
value = sprintf('%.2f ± %.2f', mean(values), std(values));
end


function value = median_iqr(values)
values = numeric_vector(values);
values = values(isfinite(values));
assert(~isempty(values), 'Cannot summarize an empty clinical field.');
quartiles = quantile(values, [0.25 0.75]);
value = sprintf('%.1f [%.1f-%.1f]', median(values), quartiles(1), quartiles(2));
end


function value = therapy_summary(values, line, denominator)
values = numeric_vector(values);
n = sum(values == line);
value = sprintf('%d (%.1f%%)', n, 100*n/denominator);
end


function T = add_msfc_sdmt(T)
required = {'SubjectID','Group','T25FW','x9HPTD','x9HPTND','SDMTcorrect'};
assert(isempty(setdiff(required, T.Properties.VariableNames)), ...
    'MSFC source columns are missing.');
subjects = string(T.SubjectID);
assert(numel(unique(subjects)) == height(T), ...
    'MSFC scoring requires one row per unique subject.');

components = [numeric_vector(T.T25FW), numeric_vector(T.x9HPTD), ...
    numeric_vector(T.x9HPTND), numeric_vector(T.SDMTcorrect)];
valid = all(isfinite(components),2) & all(components(:,1:3) > 0,2) & ...
    components(:,4) >= 0;
reference = valid & ~startsWith(subjects, 'sub-C');
assert(height(T) == 132 && sum(valid) == 127 && sum(reference) == 84, ...
    'Unexpected MSFC-SDMT cohort.');

arm = nan(height(T),1);
arm(valid) = (1./components(valid,2) + 1./components(valid,3))/2;
z_arm = nan(height(T),1);
z_leg = nan(height(T),1);
z_cog = nan(height(T),1);
z_arm(valid) = (arm(valid)-mean(arm(reference)))/std(arm(reference));
z_leg(valid) = -(components(valid,1)-mean(components(reference,1))) / ...
    std(components(reference,1));
z_cog(valid) = (components(valid,4)-mean(components(reference,4))) / ...
    std(components(reference,4));
score = mean([z_arm,z_leg,z_cog],2,'omitmissing');
score(~valid) = NaN;
T.MSFC_SDMT = score;
end


function values = numeric_vector(values)
if isnumeric(values) || islogical(values)
    values = double(values);
elseif iscell(values)
    values = str2double(string(values));
else
    values = str2double(string(values));
end
values = values(:);
end


function value = scalar_number(value)
value = numeric_vector(value);
assert(numel(value) == 1 && isfinite(value), 'Expected one finite number.');
value = value(1);
end


function write_plain_csv(path, values)
lines = strings(size(values,1),1);
for row = 1:size(values,1)
    fields = cell(1,size(values,2));
    for column = 1:size(values,2)
        raw = char(string(values{row,column}));
        needs_quotes = contains(raw, ',') || contains(raw, '"') || ...
            contains(raw, newline) || contains(raw, char(13));
        raw = strrep(raw, '"', '""');
        if needs_quotes
            raw = ['"' raw '"']; %#ok<AGROW>
        end
        fields{column} = raw;
    end
    lines(row) = string(strjoin(fields, ','));
end
writelines(lines, path, 'Encoding', 'UTF-8');
end


function path = require_file(path, label)
assert(isfile(path), '%s not found: %s', label, path);
end


function result = is_text_scalar(value)
result = ischar(value) || (isstring(value) && isscalar(value));
end
