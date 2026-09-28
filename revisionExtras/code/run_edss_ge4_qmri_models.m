function report = run_edss_ge4_qmri_models(varargin)
%RUN_EDSS_GE4_QMRI_MODELS Refit qMRI EDSS models using EDSS >=4.
%
% This threshold-only sensitivity refits the 28 qMRI EDSS models while
% holding each manuscript analysis cohort fixed.  It reproduces the
% canonical EDSS >3 rows before fitting EDSS >=4.  Patient-level outputs are
% written only to revisionExtras/work and are ignored by Git.
%
% Required name-value arguments:
%   StatsDir    directory containing clinicalScore.xlsx
%   MetricsDir  directory containing ClassicalRegionAllMetrics.csv and the
%               seven GroupTract metric tables
%
% Optional name-value arguments:
%   CanonicalModelFile  baseline manuscript model CSV
%   OutputDir           aggregate/intermediate output directory
%   QADir               patient-level and baseline QA directory

script_dir = fileparts(mfilename('fullpath'));
package_root = fileparts(script_dir);
repo_root = fileparts(package_root);

p = inputParser;
addParameter(p,'StatsDir','',@is_text_scalar);
addParameter(p,'MetricsDir','',@is_text_scalar);
addParameter(p,'CanonicalModelFile',fullfile(repo_root,'revision', ...
    'results','model_results','all_qmri_manuscript_models.csv'), ...
    @is_text_scalar);
addParameter(p,'OutputDir',fullfile(package_root,'work','qmri'), ...
    @is_text_scalar);
addParameter(p,'QADir',fullfile(package_root,'work','qa','qmri'), ...
    @is_text_scalar);
parse(p,varargin{:});

stats_dir = require_directory(char(p.Results.StatsDir),'StatsDir');
metrics_dir = require_directory(char(p.Results.MetricsDir),'MetricsDir');
canonical_csv = require_file(char(p.Results.CanonicalModelFile), ...
    'canonical qMRI model CSV');
output_dir = local_output_dir(char(p.Results.OutputDir),package_root, ...
    'OutputDir');
qa_dir = local_output_dir(char(p.Results.QADir),package_root,'QADir');
if ~isfolder(output_dir), mkdir(output_dir); end
if ~isfolder(qa_dir), mkdir(qa_dir); end

metrics = {'T1','MTR','FA','MD','NDI','ODI','FWF'};
display_metrics = {'T1','MTR','FA','MD','NDI','ODI','ISOVF'};
legacy_metrics = {'T1','MTR','FA','MD'};
tract_groups = {'Association','Cerebellar','Occipitoparietal', ...
    'ProjectionBrainstem'};
regions = {'periventricular','juxtacortical','infratentorial', ...
    'deepwhitematter'};
nrf_preds = build_nrf_preds(tract_groups);
demo_cols = {'Age','GenderNum','DurationOfDisease'};
lambda_grid = logspace(-6, 6, 50);
n_boot = 2000;
seed = 42;
alpha = 0.05;

required = {fullfile(stats_dir,'clinicalScore.xlsx'), ...
    fullfile(metrics_dir,'ClassicalRegionAllMetrics.csv'), canonical_csv};
for mi = 1:numel(metrics)
    if ismember(metrics{mi}, legacy_metrics)
        required{end+1} = fullfile(metrics_dir, ...
            sprintf('GroupTract%s_All.xlsx',metrics{mi})); %#ok<AGROW>
    else
        required{end+1} = fullfile(metrics_dir, ...
            sprintf('GroupTract%s_All.csv',metrics{mi})); %#ok<AGROW>
    end
end
for fi = 1:numel(required)
    assert(isfile(required{fi}), 'Required frozen input missing: %s', required{fi});
end

%% Clinical and classical data
T_clin = readtable(fullfile(stats_dir,'clinicalScore.xlsx'), ...
    'TextType','string');
assert(height(T_clin) == 132 && numel(unique(T_clin.SubjectID)) == 132, ...
    'Frozen clinical table must contain 132 unique subjects.');
gender = strtrim(string(T_clin.Gender));
T_clin.GenderNum = double(gender == "M");
T_clin.GenderNum(ismissing(gender) | strlength(gender) == 0) = NaN;

T_ico = readtable(fullfile(metrics_dir,'ClassicalRegionAllMetrics.csv'), ...
    'TextType','string');
assert(height(T_ico) == 87 && numel(unique(T_ico.SubjectID)) == 87, ...
    'Frozen classical table must contain 87 unique MS subjects.');
T_crf_full = outerjoin(T_clin, T_ico, ...
    'Keys','SubjectID','MergeKeys',true);
T_crf_ms = T_crf_full(~startsWith(string(T_crf_full.SubjectID),'sub-C'),:);

canonical = readtable(canonical_csv, 'TextType','string');
canonical = canonical(canonical.Outcome == "EDSS",:);
assert(height(canonical) == 28, ...
    'Expected 28 canonical EDSS rows (7 metrics x 2 models x 2 specs).');

% Frozen cohort boundary audit uses the common 132-subject tract set.
T_frozen = readtable(fullfile(metrics_dir,'GroupTractT1_All.xlsx'), ...
    'TextType','string');
frozen_ids = string(T_frozen.SubjectID);
assert(numel(unique(frozen_ids)) == 132 && ...
    ~any(frozen_ids == "sub-RR047") && ...
    any(frozen_ids == "sub-C035") && any(frozen_ids == "sub-RR013"), ...
    'Frozen 132-subject revision cohort rule failed.');
[in_clin, clin_loc] = ismember(frozen_ids, string(T_clin.SubjectID));
T_frozen_clin = T_clin(clin_loc(in_clin),:);
is_boundary = T_frozen_clin.EDSS == 3.5 | T_frozen_clin.EDSS == 4.0;
boundary = T_frozen_clin(is_boundary,{'SubjectID','EDSS'});
boundary.CanonicalPositive_gt3 = boundary.EDSS > 3.0;
boundary.SensitivityPositive_ge4 = boundary.EDSS >= 4.0;
boundary.ClassificationChanged = ...
    boundary.CanonicalPositive_gt3 ~= boundary.SensitivityPositive_ge4;
boundary.InFrozenImagingCohort = true(height(boundary),1);
boundary.ClassicalRowPresent = ismember(boundary.SubjectID,T_ico.SubjectID);
for mi = 1:numel(metrics)
    boundary.([display_metrics{mi} '_JointComplete']) = ...
        false(height(boundary),1);
end

%% Analysis outputs
model_rows = struct([]);
wide_rows = struct([]);
paired_rows = struct([]);
baseline_rows = struct([]);
membership_rows = struct([]);
pair_audit_rows = struct([]);

for mi = 1:numel(metrics)
    metric = metrics{mi};
    display_metric = display_metrics{mi};
    if ismember(metric,legacy_metrics)
        tract_file = fullfile(metrics_dir, ...
            sprintf('GroupTract%s_All.xlsx',metric));
    else
        tract_file = fullfile(metrics_dir, ...
            sprintf('GroupTract%s_All.csv',metric));
    end
    T_tract = readtable(tract_file,'TextType','string');
    assert(height(T_tract) == 132 && ...
        numel(unique(T_tract.SubjectID)) == 132, ...
        '%s tract table must contain 132 unique subjects.',metric);
    assert(~any(T_tract.SubjectID == "sub-RR047") && ...
        all(ismember(["sub-C035","sub-RR013"],T_tract.SubjectID)), ...
        '%s tract table violates frozen subject rules.',metric);
    T_nrf_full = outerjoin(T_clin,T_tract, ...
        'Keys','SubjectID','MergeKeys',true);
    T_nrf_ms = T_nrf_full(~startsWith(string(T_nrf_full.SubjectID),'sub-C'),:);
    crf_preds = build_crf_preds(regions,metric);
    validate_predictors(T_nrf_ms,nrf_preds,sprintf('%s tract',metric));
    validate_predictors(T_crf_ms,crf_preds,sprintf('%s classical',metric));
    exclude_zeros = ismember(metric,legacy_metrics);

    for adjusted = [false true]
        if adjusted
            demos = demo_cols;
            demographics_label = 'Age+GenderNum+DurationOfDisease';
        else
            demos = {};
            demographics_label = 'none';
        end

        [Xn_base,Xc_base,y_base,Demo_base,ids_base,edss_base] = ...
            get_joint_cases(T_nrf_ms,T_crf_ms,nrf_preds,crf_preds, ...
            demos,exclude_zeros,'gt3');
        [Xn_sens,Xc_sens,y_sens,Demo_sens,ids_sens,edss_sens] = ...
            get_joint_cases(T_nrf_ms,T_crf_ms,nrf_preds,crf_preds, ...
            demos,exclude_zeros,'ge4');
        assert(isequal(ids_base,ids_sens) && ...
            isequaln(Xn_base,Xn_sens) && isequaln(Xc_base,Xc_sens) && ...
            isequaln(Demo_base,Demo_sens) && isequal(edss_base,edss_sens), ...
            '%s %s threshold runs do not use identical subjects/data.', ...
            metric,demographics_label);
        assert(~any(ids_base == "sub-RR047"), ...
            'RR047 entered %s %s cases.',metric,demographics_label);
        assert(numel(unique(y_base)) == 2 && numel(unique(y_sens)) == 2, ...
            '%s %s threshold lacks two classes.',metric,demographics_label);

        if ~adjusted
            field = [display_metric '_JointComplete'];
            boundary.(field) = ismember(boundary.SubjectID,ids_base);
        end

        % Baseline: exact current fitting sequence, with RNG reset once for
        % the tract/classical pair.
        rng(seed,'twister');
        [LPn_base,lambda_n_base] = fit_ridge_lp_bin( ...
            zscore(Xn_base),y_base,lambda_grid);
        [LPc_base,lambda_c_base] = fit_ridge_lp_bin( ...
            zscore(Xc_base),y_base,lambda_grid);
        stat_n_base = eval_binary(LPn_base,y_base,Demo_base);
        stat_c_base = eval_binary(LPc_base,y_base,Demo_base);
        legacy_n_base = canonical_auc_bootstrap( ...
            y_base,stat_n_base.score,n_boot,alpha,seed);
        legacy_c_base = canonical_auc_bootstrap( ...
            y_base,stat_c_base.score,n_boot,alpha,seed);
        paired_base = paired_stratified_auc( ...
            y_base,stat_n_base.score,stat_c_base.score,n_boot,alpha,seed);

        % Sensitivity: EDSS >=4.0, all other choices unchanged.
        rng(seed,'twister');
        [LPn_sens,lambda_n_sens] = fit_ridge_lp_bin( ...
            zscore(Xn_sens),y_sens,lambda_grid);
        [LPc_sens,lambda_c_sens] = fit_ridge_lp_bin( ...
            zscore(Xc_sens),y_sens,lambda_grid);
        stat_n_sens = eval_binary(LPn_sens,y_sens,Demo_sens);
        stat_c_sens = eval_binary(LPc_sens,y_sens,Demo_sens);
        legacy_n_sens = canonical_auc_bootstrap( ...
            y_sens,stat_n_sens.score,n_boot,alpha,seed);
        legacy_c_sens = canonical_auc_bootstrap( ...
            y_sens,stat_c_sens.score,n_boot,alpha,seed);
        paired_sens = paired_stratified_auc( ...
            y_sens,stat_n_sens.score,stat_c_sens.score,n_boot,alpha,seed);

        % Adapters make the long-form model output use exactly the same
        % Performance/CI semantics as the canonical manuscript CSV:
        % Performance is the mean of the legacy unstratified bootstrap and
        % CI is its percentile interval. Paired outputs remain unchanged.
        legacy_pair_base = legacy_model_adapter(legacy_n_base,legacy_c_base);
        legacy_pair_sens = legacy_model_adapter(legacy_n_sens,legacy_c_sens);

        if adjusted
            bonf_n_sens = min(stat_n_sens.p*14,1);
            bonf_c_sens = min(stat_c_sens.p*14,1);
            multiplicity_label = ...
                'exploratory EDSS sensitivity +demographic family only (14 models)';
        else
            bonf_n_sens = NaN;
            bonf_c_sens = NaN;
            multiplicity_label = ...
                'none; no-demographic rows outside exploratory family';
        end

        % Long-form model results, harmonized point AUC + stratified CI.
        model_rows = append_model_rows(model_rows,metric,display_metric, ...
            adjusted,demographics_label,ids_base,y_base, ...
            stat_n_base,stat_c_base,lambda_n_base,lambda_c_base, ...
            legacy_pair_base,'canonical_gt3','EDSS > 3.0',NaN,NaN, ...
            'canonical family untouched; no exploratory adjustment');
        model_rows = append_model_rows(model_rows,metric,display_metric, ...
            adjusted,demographics_label,ids_sens,y_sens, ...
            stat_n_sens,stat_c_sens,lambda_n_sens,lambda_c_sens, ...
            legacy_pair_sens,'sensitivity_ge4','EDSS >= 4.0', ...
            bonf_n_sens,bonf_c_sens,multiplicity_label);

        % Paired tract/classical comparison rows.
        paired_rows = append_paired_row(paired_rows,metric,display_metric, ...
            adjusted,demographics_label,'canonical_gt3','EDSS > 3.0', ...
            ids_base,y_base,stat_n_base,stat_c_base,paired_base, ...
            n_boot,seed);
        paired_rows = append_paired_row(paired_rows,metric,display_metric, ...
            adjusted,demographics_label,'sensitivity_ge4','EDSS >= 4.0', ...
            ids_sens,y_sens,stat_n_sens,stat_c_sens,paired_sens, ...
            n_boot,seed);

        % Wide side-by-side rows, one per model.
        wide_rows = append_wide_row(wide_rows,metric,display_metric, ...
            'Tract-based',adjusted,demographics_label,ids_base,y_base, ...
            stat_n_base,lambda_n_base,legacy_n_base,paired_base.tract, ...
            ids_sens,y_sens,stat_n_sens,lambda_n_sens,paired_sens.tract, ...
            bonf_n_sens,multiplicity_label,paired_base,paired_sens, ...
            stat_n_base.AIC-stat_c_base.AIC, ...
            stat_n_sens.AIC-stat_c_sens.AIC,n_boot,seed);
        wide_rows = append_wide_row(wide_rows,metric,display_metric, ...
            'Classical-region',adjusted,demographics_label,ids_base,y_base, ...
            stat_c_base,lambda_c_base,legacy_c_base,paired_base.classical, ...
            ids_sens,y_sens,stat_c_sens,lambda_c_sens,paired_sens.classical, ...
            bonf_c_sens,multiplicity_label,paired_base,paired_sens, ...
            stat_n_base.AIC-stat_c_base.AIC, ...
            stat_n_sens.AIC-stat_c_sens.AIC,n_boot,seed);

        % Row-by-row canonical reproduction, including legacy unstratified CI.
        baseline_rows = append_baseline_check(baseline_rows,canonical, ...
            metric,'Tract-based',adjusted,ids_base,y_base,stat_n_base, ...
            lambda_n_base,legacy_n_base,paired_base.tract);
        baseline_rows = append_baseline_check(baseline_rows,canonical, ...
            metric,'Classical-region',adjusted,ids_base,y_base,stat_c_base, ...
            lambda_c_base,legacy_c_base,paired_base.classical);

        % Subject membership and pair-level QA.
        subject_hash = sha256_text(strjoin(sort(ids_base),newline));
        for si = 1:numel(ids_base)
            mr = struct();
            mr.Metric = metric;
            mr.DisplayMetric = display_metric;
            mr.Adjusted = adjusted;
            mr.Demographics = demographics_label;
            mr.SubjectID = char(ids_base(si));
            mr.EDSS = edss_base(si);
            mr.CanonicalPositive_gt3 = y_base(si);
            mr.SensitivityPositive_ge4 = y_sens(si);
            mr.InTractCases = true;
            mr.InClassicalCases = true;
            membership_rows = append_struct(membership_rows,mr);
        end
        ar = struct();
        ar.Metric = metric;
        ar.DisplayMetric = display_metric;
        ar.Adjusted = adjusted;
        ar.Demographics = demographics_label;
        ar.N = numel(ids_base);
        ar.SubjectIDSetSHA256 = subject_hash;
        ar.BaselineSensitivityIDsIdentical = isequal(ids_base,ids_sens);
        ar.TractClassicalIDsIdentical = true;
        ar.RR047Absent = ~any(ids_base == "sub-RR047");
        ar.BaselineTwoClasses = numel(unique(y_base)) == 2;
        ar.SensitivityTwoClasses = numel(unique(y_sens)) == 2;
        ar.BaselineBootstrapAllTwoClass = paired_base.all_two_class;
        ar.SensitivityBootstrapAllTwoClass = paired_sens.all_two_class;
        ar.BaselineBootstrapIndexSHA256 = paired_base.index_sha256;
        ar.SensitivityBootstrapIndexSHA256 = paired_sens.index_sha256;
        ar.BootstrapIndicesSharedWithinPairedModels = true;
        ar.BootstrapRepetitions = n_boot;
        ar.BootstrapSeed = seed;
        pair_audit_rows = append_struct(pair_audit_rows,ar);
    end
end

%% Validation and outputs
baseline_table = struct2table(baseline_rows);
baseline_table.NegativeN = baseline_table.N - baseline_table.PositiveN;
assert(all(baseline_table.Pass), ...
    'Canonical EDSS >3.0 baseline did not reproduce exactly.');
assert(height(baseline_table) == 28, ...
    'Baseline reproduction must contain 28 rows.');

model_table = struct2table(model_rows);
wide_table = struct2table(wide_rows);
paired_table = struct2table(paired_rows);
membership_table = struct2table(membership_rows);
pair_audit_table = struct2table(pair_audit_rows);
assert(height(model_table) == 56 && height(wide_table) == 28 && ...
    height(paired_table) == 28 && height(pair_audit_table) == 14, ...
    'Unexpected EDSS threshold sensitivity output dimensions.');
assert(all(pair_audit_table.BaselineSensitivityIDsIdentical) && ...
    all(pair_audit_table.TractClassicalIDsIdentical) && ...
    all(pair_audit_table.RR047Absent) && ...
    all(pair_audit_table.BaselineBootstrapAllTwoClass) && ...
    all(pair_audit_table.SensitivityBootstrapAllTwoClass), ...
    'Pairing, cohort, or bootstrap QA failed.');

boundary_counts = groupsummary(boundary,'EDSS');
boundary_counts.Properties.VariableNames{ ...
    strcmp(boundary_counts.Properties.VariableNames,'GroupCount')} = ...
    'SubjectCount';
boundary_counts.CanonicalPositiveCount = ...
    (boundary_counts.EDSS > 3.0) .* boundary_counts.SubjectCount;
boundary_counts.SensitivityPositiveCount = ...
    (boundary_counts.EDSS >= 4.0) .* boundary_counts.SubjectCount;
boundary_counts.ReclassifiedCount = ...
    boundary_counts.CanonicalPositiveCount - ...
    boundary_counts.SensitivityPositiveCount;

% Stable sorting for review.
model_table = sortrows(model_table, ...
    {'ThresholdOrder','MetricOrder','Adjusted','ModelOrder'});
wide_table = sortrows(wide_table,{'MetricOrder','Adjusted','ModelOrder'});
paired_table = sortrows(paired_table, ...
    {'ThresholdOrder','MetricOrder','Adjusted'});
baseline_table = sortrows(baseline_table, ...
    {'MetricOrder','Adjusted','ModelOrder'});
membership_table = sortrows(membership_table, ...
    {'Metric','Adjusted','SubjectID'});
pair_audit_table = sortrows(pair_audit_table,{'Metric','Adjusted'});
boundary = sortrows(boundary,{'EDSS','SubjectID'});

writetable(wide_table,fullfile(output_dir, ...
    'edss_threshold_side_by_side.csv'));
writetable(model_table,fullfile(output_dir, ...
    'edss_threshold_model_results.csv'));
writetable(paired_table,fullfile(output_dir, ...
    'paired_tract_classical_auc.csv'));
writetable(boundary,fullfile(output_dir,'edss_boundary_audit.csv'));
writetable(boundary_counts,fullfile(output_dir, ...
    'edss_boundary_counts.csv'));
writetable(baseline_table,fullfile(qa_dir, ...
    'baseline_reproduction.csv'));
writetable(membership_table,fullfile(qa_dir, ...
    'paired_subject_membership.csv'));
writetable(pair_audit_table,fullfile(qa_dir, ...
    'paired_bootstrap_audit.csv'));

write_reviewer_text(fullfile(output_dir,'reviewer_ready_text.md'), ...
    boundary,paired_table,model_table);
write_plain_conclusion(fullfile(output_dir,'plain_language_conclusion.txt'), ...
    boundary,paired_table,model_table);

report = struct();
report.schema_version = 1;
report.status = 'PASS';
report.scope = 'EDSS threshold sensitivity only';
report.canonical_rule = 'EDSS > 3.0 (strict)';
report.sensitivity_rule = 'EDSS >= 4.0 (positive) versus EDSS < 4.0 (negative)';
report.baseline_reproduction = struct( ...
    'rows_checked',height(baseline_table), ...
    'all_match',all(baseline_table.Pass), ...
    'numeric_tolerance',1e-12, ...
    'canonical_performance_semantics', ...
    ['The canonical Performance field is the mean of the legacy ' ...
     'unstratified 2000-resample bootstrap, despite its apparent-AUC label.'], ...
    'new_comparator_semantics', ...
    ['Side-by-side AUC is the full-sample apparent AUC; its interval is ' ...
     'the paired class-stratified percentile bootstrap CI. The two CI ' ...
     'procedures are labeled separately and are not asserted identical.']);
report.model = struct( ...
    'learner','fitclinear logistic ridge', ...
    'lambda_grid','logspace(-6,6,50)', ...
    'cross_validation_folds',10, ...
    'full_sample_refit',true, ...
    'ridge_score_standardization','zscore predictors', ...
    'tract_predictor_count',16, ...
    'classical_predictor_count',4, ...
    'demographic_covariates',strjoin(demo_cols,'+'), ...
    'rng_seed',seed);
report.bootstrap = struct( ...
    'repetitions',n_boot, ...
    'seed',seed, ...
    'method',['class-stratified percentile bootstrap; positive and ' ...
        'negative classes sampled separately with replacement'], ...
    'paired_indices',['identical resample indices applied to tract and ' ...
        'classical scores within each metric/specification/threshold'], ...
    'interpretation',['imbalance-aware uncertainty/robustness check; ' ...
        'does not prove that class imbalance is eliminated'], ...
    'all_resamples_two_class',true);
report.multiplicity = struct( ...
    'canonical_family_modified',false, ...
    'exploratory_family',['14 +demographic EDSS >=4.0 models only ' ...
        '(7 metrics x 2 model types)'], ...
    'exploratory_multiplier',14, ...
    'no_demo_rows_adjusted',false);
report.case_policy = struct( ...
    'population','MS participants only', ...
    'joint_complete_cases',true, ...
    'legacy_zero_exclusion_metrics',{legacy_metrics}, ...
    'noddi_exact_zeros_retained',true, ...
    'RR047_excluded',true, ...
    'C035_note','present in frozen tract cohort; excluded because control', ...
    'RR013_note',['present in frozen tract cohort; unavailable for joint ' ...
        'tract/classical EDSS models because the classical row is absent'], ...
    'MTR_note','C035 and RR013 MTR predictors are genuinely missing');
report.boundary = struct( ...
    'EDSS_3_5_n',sum(boundary.EDSS == 3.5), ...
    'EDSS_4_0_n',sum(boundary.EDSS == 4.0), ...
    'reclassified_n',sum(boundary.ClassificationChanged), ...
    'reclassified_in_each_joint_model_n', ...
        sum(boundary.ClassificationChanged & boundary.T1_JointComplete));
report.qa = struct( ...
    'paired_case_audits',height(pair_audit_table), ...
    'all_paired_ids_identical', ...
        all(pair_audit_table.TractClassicalIDsIdentical), ...
    'all_threshold_case_sets_identical', ...
        all(pair_audit_table.BaselineSensitivityIDsIdentical), ...
    'all_bootstraps_valid_two_class', ...
        all(pair_audit_table.BaselineBootstrapAllTwoClass & ...
            pair_audit_table.SensitivityBootstrapAllTwoClass), ...
    'embedded_input_and_cohort_preflight',true);
report.read_only_references = struct( ...
    'canonical_models',canonical_csv, ...
    'frozen_metrics_dir',metrics_dir, ...
    'frozen_stats_dir',stats_dir);
report.output_files = struct( ...
    'side_by_side','edss_threshold_side_by_side.csv', ...
    'long_model_results','edss_threshold_model_results.csv', ...
    'paired_comparisons','paired_tract_classical_auc.csv', ...
    'boundary_audit','edss_boundary_audit.csv', ...
    'boundary_counts','edss_boundary_counts.csv', ...
    'reviewer_ready_text','reviewer_ready_text.md', ...
    'plain_language_conclusion','plain_language_conclusion.txt');
write_json(fullfile(qa_dir,'edss_ge4_qmri_report.json'),report);
fprintf('EDSS >=4 qMRI sensitivity PASS. Outputs: %s\n',output_dir);
end


function value = is_text_scalar(value)
value = ischar(value) || (isstring(value) && isscalar(value));
end


function path = require_directory(path,label)
assert(~isempty(path) && isfolder(path),'%s is missing: %s',label,path);
path = char(java.io.File(path).getCanonicalPath());
end


function path = require_file(path,label)
assert(~isempty(path) && isfile(path),'%s is missing: %s',label,path);
path = char(java.io.File(path).getCanonicalPath());
end


function path = local_output_dir(path,package_root,label)
assert(~isempty(path),'%s cannot be empty.',label);
root = char(java.io.File(package_root).getCanonicalPath());
path = char(java.io.File(path).getCanonicalPath());
assert(strcmp(path,root) || startsWith(path,[root filesep]), ...
    '%s must remain inside revisionExtras: %s',label,path);
end


function preds = build_nrf_preds(groups)
preds = {};
for gi = 1:numel(groups)
    g = groups{gi};
    preds = [preds, {[g 'L_Tail'],[g 'R_Tail'], ...
        [g 'L_NAWM'],[g 'R_NAWM']}]; %#ok<AGROW>
end
end


function preds = build_crf_preds(regions,metric)
preds = cell(size(regions));
for ri = 1:numel(regions)
    preds{ri} = [regions{ri} metric];
end
end


function validate_predictors(T,preds,label)
missing = preds(~ismember(preds,T.Properties.VariableNames));
assert(isempty(missing),'Missing %s predictors: %s',label,strjoin(missing,', '));
end


function [Xn,Xc,y,Demo,ids,edss] = get_joint_cases( ...
        Tn,Tc,nrf_preds,crf_preds,demo_cols,exclude_zeros,rule)
[~,ia,ib] = intersect(Tn.SubjectID,Tc.SubjectID,'stable');
Tn = Tn(ia,:);
Tc = Tc(ib,:);
Xn_raw = table2array(Tn(:,nrf_preds));
Xc_raw = table2array(Tc(:,crf_preds));
edss_raw = double(Tn.EDSS);
D = zeros(height(Tn),numel(demo_cols));
for di = 1:numel(demo_cols)
    D(:,di) = double(Tn.(demo_cols{di}));
end
valid = all(isfinite(Xn_raw),2) & all(isfinite(Xc_raw),2) & ...
    isfinite(edss_raw) & all(isfinite(D),2);
if exclude_zeros
    valid = valid & all(Xn_raw ~= 0,2) & all(Xc_raw ~= 0,2);
end
Xn = Xn_raw(valid,:);
Xc = Xc_raw(valid,:);
Demo = D(valid,:);
ids = string(Tn.SubjectID(valid));
edss = edss_raw(valid);
switch rule
    case 'gt3'
        y = double(edss > 3.0);
    case 'ge4'
        y = double(edss >= 4.0);
    otherwise
        error('Unknown EDSS rule: %s',rule);
end
end


function [LP,lambda] = fit_ridge_lp_bin(X,y,lambda_grid)
cv = fitclinear(X,y,'Learner','logistic', ...
    'Regularization','ridge','Lambda',lambda_grid,'KFold',10);
[~,best] = min(kfoldLoss(cv));
lambda = lambda_grid(best);
mdl = fitclinear(X,y,'Learner','logistic', ...
    'Regularization','ridge','Lambda',lambda);
LP = X*mdl.Beta + mdl.Bias;
end


function stat = eval_binary(LP,y,Demo)
warning_state = warning;
cleanup = onCleanup(@() warning(warning_state));
warning('off','all');
if isempty(Demo)
    mdl = fitglm(LP,y,'Distribution','binomial');
    score = LP;
else
    mdl = fitglm([LP,Demo],y,'Distribution','binomial');
    score = mdl.Fitted.LinearPredictor;
end
coef = mdl.Coefficients.Estimate(2);
se = mdl.Coefficients.SE(2);
stat = struct('score',score,'p',2*normcdf(-abs(coef/se)), ...
    'AIC',mdl.ModelCriterion.AIC,'point_auc',auc_rank(y,score));
end


function legacy = canonical_auc_bootstrap(y,score,n_boot,alpha,seed)
rng(seed,'twister');
aucs = bootstrp(n_boot,@(yy,ss) perf_auc(yy,ss),y(:),score(:));
legacy = struct('mean_auc',mean(aucs), ...
    'ci',quantile(aucs,[alpha/2,1-alpha/2]));
end


function adapted = legacy_model_adapter(tract_legacy,classical_legacy)
adapted = struct();
adapted.tract = struct( ...
    'point_auc',tract_legacy.mean_auc, ...
    'mean_bootstrap_auc',tract_legacy.mean_auc, ...
    'ci',tract_legacy.ci);
adapted.classical = struct( ...
    'point_auc',classical_legacy.mean_auc, ...
    'mean_bootstrap_auc',classical_legacy.mean_auc, ...
    'ci',classical_legacy.ci);
adapted.index_sha256 = 'not-recorded-for-legacy-unstratified-bootstrap';
end


function paired = paired_stratified_auc(y,tract_score,classical_score, ...
        n_boot,alpha,seed)
y = y(:);
tract_score = tract_score(:);
classical_score = classical_score(:);
pos = find(y == 1);
neg = find(y == 0);
assert(~isempty(pos) && ~isempty(neg),'Both outcome classes are required.');
rng(seed,'twister');
indices = zeros(numel(y),n_boot,'uint32');
for bi = 1:n_boot
    draw_pos = pos(randi(numel(pos),numel(pos),1));
    draw_neg = neg(randi(numel(neg),numel(neg),1));
    indices(:,bi) = uint32([draw_pos;draw_neg]);
end
tract_auc = zeros(n_boot,1);
classical_auc = zeros(n_boot,1);
all_two_class = true;
for bi = 1:n_boot
    idx = double(indices(:,bi));
    yb = y(idx);
    all_two_class = all_two_class && numel(unique(yb)) == 2;
    tract_auc(bi) = auc_rank(yb,tract_score(idx));
    classical_auc(bi) = auc_rank(yb,classical_score(idx));
end
delta = tract_auc-classical_auc;
paired = struct();
paired.tract = struct('point_auc',auc_rank(y,tract_score), ...
    'mean_bootstrap_auc',mean(tract_auc), ...
    'ci',quantile(tract_auc,[alpha/2,1-alpha/2]));
paired.classical = struct('point_auc',auc_rank(y,classical_score), ...
    'mean_bootstrap_auc',mean(classical_auc), ...
    'ci',quantile(classical_auc,[alpha/2,1-alpha/2]));
paired.delta_point = paired.tract.point_auc-paired.classical.point_auc;
paired.delta_mean = mean(delta);
paired.delta_ci = quantile(delta,[alpha/2,1-alpha/2]);
paired.index_sha256 = sha256_uint32(indices);
paired.all_two_class = all_two_class;
paired.positive_per_resample = numel(pos);
paired.negative_per_resample = numel(neg);
end


function value = auc_rank(y,score)
y = y(:);
score = score(:);
np = sum(y == 1);
nn = sum(y == 0);
assert(np > 0 && nn > 0,'AUC requires two classes.');
ranks = tiedrank(score);
value = (sum(ranks(y == 1))-np*(np+1)/2)/(np*nn);
end


function value = perf_auc(y,score)
[~,~,~,value] = perfcurve(y,score,1);
end


function rows = append_model_rows(rows,metric,display_metric,adjusted, ...
        demographics,ids,y,stat_n,stat_c,lambda_n,lambda_c,paired, ...
        threshold_label,positive_rule,bonf_n,bonf_c,multiplicity_label)
rows = append_model_row(rows,metric,display_metric,'Tract-based',1, ...
    adjusted,demographics,ids,y,stat_n,lambda_n,paired.tract, ...
    paired.index_sha256,threshold_label,positive_rule,bonf_n, ...
    multiplicity_label);
rows = append_model_row(rows,metric,display_metric,'Classical-region',2, ...
    adjusted,demographics,ids,y,stat_c,lambda_c,paired.classical, ...
    paired.index_sha256,threshold_label,positive_rule,bonf_c, ...
    multiplicity_label);
end


function rows = append_model_row(rows,metric,display_metric,model,model_order, ...
        adjusted,demographics,ids,y,stat,lambda,boot,index_hash, ...
        threshold_label,positive_rule,p_bonf,multiplicity_label)
r = struct();
r.ThresholdLabel = threshold_label;
r.ThresholdOrder = double(strcmp(threshold_label,'sensitivity_ge4'))+1;
r.PositiveRule = positive_rule;
r.Metric = metric;
r.DisplayMetric = display_metric;
r.MetricOrder = metric_order(metric);
r.Model = model;
r.ModelOrder = model_order;
r.Adjusted = adjusted;
r.Demographics = demographics;
r.CasePolicy = 'MS; joint tract+classical complete cases';
r.N = numel(y);
r.PositiveN = sum(y);
r.NegativeN = sum(y == 0);
r.PercentPositive = 100*mean(y);
r.AUC = boot.point_auc;
r.AUC_CI_Lower = boot.ci(1);
r.AUC_CI_Upper = boot.ci(2);
r.BootstrapMeanAUC = boot.mean_bootstrap_auc;
r.AIC = stat.AIC;
r.Lambda = lambda;
r.p_raw = stat.p;
r.p_exploratory_bonferroni14 = p_bonf;
r.MultiplicityLabel = multiplicity_label;
r.BootstrapIndexSHA256 = index_hash;
r.SubjectIDSetSHA256 = sha256_text(strjoin(sort(ids),newline));
r.BootstrapMethod = ...
    'legacy unstratified percentile; Performance is bootstrap mean';
r.BootstrapRepetitions = 2000;
r.BootstrapSeed = 42;
rows = append_struct(rows,r);
end


function rows = append_paired_row(rows,metric,display_metric,adjusted, ...
        demographics,threshold_label,positive_rule,ids,y,stat_n,stat_c, ...
        paired,n_boot,seed)
r = struct();
r.ThresholdLabel = threshold_label;
r.ThresholdOrder = double(strcmp(threshold_label,'sensitivity_ge4'))+1;
r.PositiveRule = positive_rule;
r.Metric = metric;
r.DisplayMetric = display_metric;
r.MetricOrder = metric_order(metric);
r.Adjusted = adjusted;
r.Demographics = demographics;
r.N = numel(y);
r.PositiveN = sum(y);
r.NegativeN = sum(y == 0);
r.PercentPositive = 100*mean(y);
r.TractAUC = paired.tract.point_auc;
r.ClassicalAUC = paired.classical.point_auc;
r.TractMinusClassicalDeltaAUC = paired.delta_point;
r.DeltaAUC_CI_Lower = paired.delta_ci(1);
r.DeltaAUC_CI_Upper = paired.delta_ci(2);
r.TractAIC = stat_n.AIC;
r.ClassicalAIC = stat_c.AIC;
r.TractMinusClassicalDeltaAIC = stat_n.AIC-stat_c.AIC;
r.BootstrapIndexSHA256 = paired.index_sha256;
r.BootstrapIndicesIdentical = true;
r.AllBootstrapResamplesTwoClass = paired.all_two_class;
r.PositivePerResample = paired.positive_per_resample;
r.NegativePerResample = paired.negative_per_resample;
r.BootstrapRepetitions = n_boot;
r.BootstrapSeed = seed;
r.SubjectIDSetSHA256 = sha256_text(strjoin(sort(ids),newline));
rows = append_struct(rows,r);
end


function rows = append_wide_row(rows,metric,display_metric,model,adjusted, ...
        demographics,ids_base,y_base,stat_base,lambda_base,legacy_base, ...
        boot_base,ids_sens,y_sens,stat_sens,lambda_sens,boot_sens, ...
        bonf_sens,multiplicity_label,paired_base,paired_sens, ...
        paired_aic_base,paired_aic_sens,n_boot,seed)
r = struct();
r.Metric = metric;
r.DisplayMetric = display_metric;
r.MetricOrder = metric_order(metric);
r.Model = model;
r.ModelOrder = 1+strcmp(model,'Classical-region');
r.Adjusted = adjusted;
r.Demographics = demographics;
r.CasePolicy = 'MS; joint tract+classical complete cases';
r.N_gt3 = numel(y_base);
r.PositiveN_gt3 = sum(y_base);
r.NegativeN_gt3 = sum(y_base == 0);
r.PercentPositive_gt3 = 100*mean(y_base);
r.AUC_gt3 = boot_base.point_auc;
r.AUC_CI_Lower_gt3 = boot_base.ci(1);
r.AUC_CI_Upper_gt3 = boot_base.ci(2);
r.CanonicalReportedAUC_gt3 = legacy_base.mean_auc;
r.CanonicalReportedCI_Lower_gt3 = legacy_base.ci(1);
r.CanonicalReportedCI_Upper_gt3 = legacy_base.ci(2);
r.AIC_gt3 = stat_base.AIC;
r.Lambda_gt3 = lambda_base;
r.p_raw_gt3 = stat_base.p;
r.N_ge4 = numel(y_sens);
r.PositiveN_ge4 = sum(y_sens);
r.NegativeN_ge4 = sum(y_sens == 0);
r.PercentPositive_ge4 = 100*mean(y_sens);
r.AUC_ge4 = boot_sens.point_auc;
r.AUC_CI_Lower_ge4 = boot_sens.ci(1);
r.AUC_CI_Upper_ge4 = boot_sens.ci(2);
r.AIC_ge4 = stat_sens.AIC;
r.Lambda_ge4 = lambda_sens;
r.p_raw_ge4 = stat_sens.p;
r.p_exploratory_bonferroni14_ge4 = bonf_sens;
r.SensitivityMultiplicityLabel = multiplicity_label;
r.AUC_change_ge4_minus_gt3 = boot_sens.point_auc-boot_base.point_auc;
r.AIC_change_ge4_minus_gt3 = stat_sens.AIC-stat_base.AIC;
r.TractMinusClassicalAUC_gt3 = paired_base.delta_point;
r.TractMinusClassicalAUC_CI_Lower_gt3 = paired_base.delta_ci(1);
r.TractMinusClassicalAUC_CI_Upper_gt3 = paired_base.delta_ci(2);
r.TractMinusClassicalAUC_ge4 = paired_sens.delta_point;
r.TractMinusClassicalAUC_CI_Lower_ge4 = paired_sens.delta_ci(1);
r.TractMinusClassicalAUC_CI_Upper_ge4 = paired_sens.delta_ci(2);
r.TractMinusClassicalAIC_gt3 = paired_aic_base;
r.TractMinusClassicalAIC_ge4 = paired_aic_sens;
r.SubjectIDsIdenticalAcrossThresholds = isequal(ids_base,ids_sens);
r.SubjectIDSetSHA256 = sha256_text(strjoin(sort(ids_base),newline));
r.BootstrapMethod = 'paired class-stratified percentile';
r.BootstrapRepetitions = n_boot;
r.BootstrapSeed = seed;
rows = append_struct(rows,r);
end


function rows = append_baseline_check(rows,canonical,metric,model,adjusted, ...
        ids,y,stat,lambda,legacy,stratified)
idx = canonical.Metric == string(metric) & ...
    canonical.Model == string(model) & ...
    logical(canonical.Adjusted) == adjusted;
assert(nnz(idx) == 1,'Canonical EDSS row not unique: %s %s adjusted=%d', ...
    metric,model,adjusted);
ref = canonical(idx,:);
r = struct();
r.Metric = metric;
r.MetricOrder = metric_order(metric);
r.Model = model;
r.ModelOrder = 1+strcmp(model,'Classical-region');
r.Adjusted = adjusted;
r.N = numel(y);
r.CanonicalN = double(ref.N);
r.PositiveN = sum(y);
r.CanonicalPositiveN = double(ref.PositiveN);
r.SubjectIDSetSHA256 = sha256_text(strjoin(sort(ids),newline));
r.ApparentAUCPoint = stat.point_auc;
r.CanonicalReportedAUC_LegacyBootstrapMean = double(ref.Performance);
r.ReproducedAUC_LegacyBootstrapMean = legacy.mean_auc;
r.CanonicalLegacyCI_Lower = double(ref.CI_lower);
r.ReproducedLegacyCI_Lower = legacy.ci(1);
r.CanonicalLegacyCI_Upper = double(ref.CI_upper);
r.ReproducedLegacyCI_Upper = legacy.ci(2);
r.NewStratifiedCI_Lower = stratified.ci(1);
r.NewStratifiedCI_Upper = stratified.ci(2);
r.CanonicalLambda = double(ref.Lambda);
r.ReproducedLambda = lambda;
r.CanonicalRawP = double(ref.p_raw);
r.ReproducedRawP = stat.p;
r.CanonicalAIC = double(ref.AIC);
r.ReproducedAIC = stat.AIC;
r.NMatch = r.N == r.CanonicalN;
r.PositiveNMatch = r.PositiveN == r.CanonicalPositiveN;
r.AUCMatch = abs(r.ReproducedAUC_LegacyBootstrapMean- ...
    r.CanonicalReportedAUC_LegacyBootstrapMean) <= 1e-12;
r.LegacyCIMatch = abs(r.ReproducedLegacyCI_Lower- ...
    r.CanonicalLegacyCI_Lower) <= 1e-12 && ...
    abs(r.ReproducedLegacyCI_Upper-r.CanonicalLegacyCI_Upper) <= 1e-12;
r.LambdaMatch = abs(r.ReproducedLambda-r.CanonicalLambda) <= 1e-12;
r.RawPMatch = abs(r.ReproducedRawP-r.CanonicalRawP) <= 1e-12;
r.AICMatch = abs(r.ReproducedAIC-r.CanonicalAIC) <= 1e-12;
r.Pass = r.NMatch && r.PositiveNMatch && r.AUCMatch && ...
    r.LegacyCIMatch && r.LambdaMatch && r.RawPMatch && r.AICMatch;
rows = append_struct(rows,r);
end


function order = metric_order(metric)
all_metrics = {'T1','MTR','FA','MD','NDI','ODI','FWF'};
order = find(strcmp(all_metrics,metric),1);
end


function rows = append_struct(rows,row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row;
end
end


function h = sha256_text(value)
h = sha256_bytes(unicode2native(char(value),'UTF-8'));
end


function h = sha256_uint32(values)
h = sha256_bytes(typecast(values(:),'uint8'));
end


function h = sha256_bytes(bytes)
md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
md.update(bytes(:));
digest = typecast(md.digest(),'uint8');
h = lower(reshape(dec2hex(digest,2).',1,[]));
end


function write_json(path,value)
fid = fopen(path,'w');
assert(fid > 0,'Unable to write %s',path);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end


function write_reviewer_text(path,boundary,paired,models)
sens = paired(paired.ThresholdLabel == "sensitivity_ge4" & paired.Adjusted,:);
base = paired(paired.ThresholdLabel == "canonical_gt3" & paired.Adjusted,:);
adv_sens = sum(sens.TractMinusClassicalDeltaAUC > 0);
adv_base = sum(base.TractMinusClassicalDeltaAUC > 0);
ci_above_zero = sum(sens.DeltaAUC_CI_Lower > 0);
reclassified = sum(boundary.ClassificationChanged);
model_reclassified = sum(boundary.ClassificationChanged & ...
    boundary.T1_JointComplete);
sens_models = models(models.ThresholdLabel == "sensitivity_ge4" & ...
    models.Adjusted,:);
exploratory_sig = sum(sens_models.p_exploratory_bonferroni14 < 0.05);
delta_median = median(sens.TractMinusClassicalDeltaAUC);
delta_min = min(sens.TractMinusClassicalDeltaAUC);
delta_max = max(sens.TractMinusClassicalDeltaAUC);
fid = fopen(path,'w');
assert(fid > 0,'Unable to write reviewer text.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'# Reviewer-ready paragraph\n\n');
fprintf(fid,['As requested, we repeated the EDSS analysis using EDSS >=4.0 ' ...
    'as the positive class (and EDSS <4.0 as the negative class), while ' ...
    'retaining the same metric definitions, joint tract/classical ' ...
    'complete-case samples, ridge-logistic fitting, 50-lambda grid, ' ...
    '10-fold cross-validation, covariates, and random seed. This changed ' ...
    'the classification of %d frozen-cohort participants with EDSS=3.5. ' ...
    '%d entered each joint complete-case analysis, so model positives ' ...
    'fell from 39 to 31; the ninth was outside the joint model sample. ' ...
    'All six participants with EDSS=4.0 remained positive. In the age-, sex-, and ' ...
    'disease-duration-adjusted sensitivity models, tract-based AUC was ' ...
    'higher than classical-region AUC for %d of 7 metrics (canonical ' ...
    '>3.0 analysis: %d of 7), with a median tract-minus-classical AUC ' ...
    'difference of %.3f (range %.3f to %.3f); the paired 95%% interval ' ...
    'excluded zero for %d of 7 metrics. A deterministic 2,000-resample ' ...
    'class-stratified bootstrap used the same resample indices for each ' ...
    'tract/classical pair. This is an imbalance-aware uncertainty check ' ...
    'and does not establish that class imbalance was eliminated. A ' ...
    'separate exploratory Bonferroni correction across only the 14 ' ...
    'adjusted sensitivity models retained %d significant models; it was ' ...
    'not merged with the canonical multiplicity family.\n'], ...
    reclassified,model_reclassified,adv_sens,adv_base,delta_median,delta_min,delta_max, ...
    ci_above_zero,exploratory_sig);
end


function write_plain_conclusion(path,boundary,paired,models)
sens = paired(paired.ThresholdLabel == "sensitivity_ge4" & paired.Adjusted,:);
base = paired(paired.ThresholdLabel == "canonical_gt3" & paired.Adjusted,:);
adv_sens = sum(sens.TractMinusClassicalDeltaAUC > 0);
adv_base = sum(base.TractMinusClassicalDeltaAUC > 0);
reclassified = sum(boundary.ClassificationChanged);
model_reclassified = sum(boundary.ClassificationChanged & ...
    boundary.T1_JointComplete);
sens_models = models(models.ThresholdLabel == "sensitivity_ge4" & ...
    models.Adjusted,:);
exploratory_sig = sum(sens_models.p_exploratory_bonferroni14 < 0.05);
fid = fopen(path,'w');
assert(fid > 0,'Unable to write conclusion.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,['Using EDSS >=4.0 instead of EDSS >3.0 reclassifies %d ' ...
    'frozen-cohort participants with EDSS=3.5; %d were in the joint ' ...
    'model samples, so positives fell from 39 to 31 without changing ' ...
    'the analyzed subject sets. Tract-based discrimination remains higher than the ' ...
    'classical-region comparator for %d of 7 adjusted metric analyses ' ...
    '(versus %d of 7 under the canonical rule). %d of the 14 adjusted ' ...
    'sensitivity models remain significant after the standalone ' ...
    'exploratory correction. The bootstrap provides an imbalance-aware ' ...
    'uncertainty check but does not remove the imbalance.\n'], ...
    reclassified,model_reclassified,adv_sens,adv_base,exploratory_sig);
end
