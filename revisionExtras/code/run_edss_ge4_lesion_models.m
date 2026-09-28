function report = run_edss_ge4_lesion_models(varargin)
%RUN_EDSS_GE4_LESION_MODELS Refit lesion-load EDSS models using EDSS >=4.
%
% Reproduces the eight canonical unadjusted lesion-load EDSS models at
% EDSS >3, then refits the same subjects and predictors at EDSS >=4.
% Reads only local inputs and writes only inside revisionExtras/work.
%
% Required name-value arguments:
%   ClinicalFile       frozen clinicalScore.xlsx
%   TractLesionFile    frozen GroupTractLesionLoad.xlsx retaining RR047
%   ClassicalLesionFile frozen icometrixLesionLoad.xlsx retaining RR047
%
% Optional name-value arguments:
%   CanonicalModelFile baseline manuscript lesion model CSV
%   OutputDir          aggregate/intermediate output directory
%   QADir              patient-level and baseline QA directory

script_dir = fileparts(mfilename('fullpath'));
package_root = fileparts(script_dir);
repo_root = fileparts(package_root);

p = inputParser;
addParameter(p,'ClinicalFile','',@is_text_scalar);
addParameter(p,'TractLesionFile','',@is_text_scalar);
addParameter(p,'ClassicalLesionFile','',@is_text_scalar);
addParameter(p,'CanonicalModelFile',fullfile(repo_root,'revision', ...
    'results','model_results','lesionload_manuscript_models.csv'), ...
    @is_text_scalar);
addParameter(p,'OutputDir',fullfile(package_root,'work','lesion'), ...
    @is_text_scalar);
addParameter(p,'QADir',fullfile(package_root,'work','qa','lesion'), ...
    @is_text_scalar);
parse(p,varargin{:});

clinical_file = require_file(char(p.Results.ClinicalFile),'ClinicalFile');
tract_file = require_file(char(p.Results.TractLesionFile), ...
    'TractLesionFile');
classical_file = require_file(char(p.Results.ClassicalLesionFile), ...
    'ClassicalLesionFile');
canonical_file = require_file(char(p.Results.CanonicalModelFile), ...
    'canonical lesion model CSV');
output_dir = local_output_dir(char(p.Results.OutputDir),package_root, ...
    'OutputDir');
qa_dir = local_output_dir(char(p.Results.QADir),package_root,'QADir');
if ~isfolder(output_dir), mkdir(output_dir); end
if ~isfolder(qa_dir), mkdir(qa_dir); end

required = {clinical_file,tract_file,classical_file,canonical_file};
for fi = 1:numel(required)
    assert(isfile(required{fi}), 'Required local input missing: %s', required{fi});
end

T_clin = readtable(clinical_file, 'TextType','string');
T_tract = readtable(tract_file, 'TextType','string');
T_ico = readtable(classical_file, 'TextType','string');
assert(height(T_clin) == 132 && numel(unique(T_clin.SubjectID)) == 132, ...
    'Clinical table must contain 132 unique subjects.');
assert(height(T_tract) == 90 && numel(unique(T_tract.SubjectID)) == 90, ...
    'Frozen tract lesion table must contain 90 unique subjects.');
assert(height(T_ico) == 90 && numel(unique(T_ico.SubjectID)) == 90, ...
    'Frozen classical lesion table must contain 90 unique subjects.');
assert(any(T_tract.SubjectID == "sub-RR047") && ...
       any(T_ico.SubjectID == "sub-RR047"), ...
    ['The frozen inputs must retain RR047 because this isolated sensitivity ' ...
     'holds the current canonical Table S5 cohort fixed.']);

T_nrf = outerjoin(T_clin,T_tract,'Keys','SubjectID','MergeKeys',true);
T_crf = outerjoin(T_clin,T_ico,'Keys','SubjectID','MergeKeys',true);
T_nrf_ms = T_nrf(~startsWith(T_nrf.SubjectID,'sub-C'),:);
T_crf_ms = T_crf(~startsWith(T_crf.SubjectID,'sub-C'),:);

tract_groups = {'Association','Cerebellar','Occipitoparietal','PB'};
regions = {'periventricular','juxtacortical','infratentorial', ...
    'deepwhitematter'};
metrics = {'LN','LV','Lnorm'};
n_boot = 2000;
alpha = 0.05;

threshold_labels = {'canonical_gt3','sensitivity_ge4'};
positive_rules = {'EDSS > 3.0','EDSS >= 4.0'};
rows = struct([]);
membership_rows = struct([]);

for ti = 1:2
    threshold_label = threshold_labels{ti};
    positive_rule = positive_rules{ti};
    use_geq = ti == 2;
    threshold = ternary(use_geq,4,3);

    % Whole-brain single-predictor models.
    wb_columns = {'WBLN','WBLV'};
    wb_metrics = {'LN','LV'};
    for wi = 1:2
        [X,y,ids,edss] = get_complete_cases(T_nrf_ms,{wb_columns{wi}}, ...
            'EDSS',threshold,use_geq);
        res = run_model_binary(X,y,n_boot,alpha);
        rows = append_result(rows,threshold_label,positive_rule, ...
            wb_metrics{wi},'Whole Brain',ids,y,res);
        membership_rows = append_membership(membership_rows, ...
            threshold_label,wb_metrics{wi},'Whole Brain',ids,edss,y);
    end

    % Joint tract/classical complete cases for LN, LV, and Lnorm.
    for mi = 1:numel(metrics)
        metric = metrics{mi};
        nrf_preds = build_nrf_preds(tract_groups,metric);
        crf_preds = build_crf_preds(regions,metric);
        validate_predictors(T_nrf_ms,nrf_preds,['Tract ' metric]);
        validate_predictors(T_crf_ms,crf_preds,['Classical ' metric]);
        [Xn,Xc,y,ids,edss] = get_joint_cases(T_nrf_ms,T_crf_ms, ...
            nrf_preds,crf_preds,'EDSS',threshold,use_geq);
        res_n = run_model_binary(Xn,y,n_boot,alpha);
        rows = append_result(rows,threshold_label,positive_rule, ...
            metric,'Tract-based',ids,y,res_n);
        membership_rows = append_membership(membership_rows, ...
            threshold_label,metric,'Tract-based',ids,edss,y);
        res_c = run_model_binary(Xc,y,n_boot,alpha);
        rows = append_result(rows,threshold_label,positive_rule, ...
            metric,'Classical-region',ids,y,res_c);
        membership_rows = append_membership(membership_rows, ...
            threshold_label,metric,'Classical-region',ids,edss,y);
    end
end

result_table = struct2table(rows);
result_table = sort_result_table(result_table);
assert(height(result_table) == 16, ...
    'Expected 16 rows (8 models x 2 thresholds).');
assert(all(result_table.N == 89), ...
    'All lesion EDSS models must retain the canonical N=89 cohort.');
assert(all(result_table.PositiveN(result_table.ThresholdLabel == ...
    "canonical_gt3") == 41), 'Canonical positive count must be 41.');
assert(all(result_table.PositiveN(result_table.ThresholdLabel == ...
    "sensitivity_ge4") == 32), 'Sensitivity positive count must be 32.');

baseline = validate_baseline(result_table,canonical_file);
assert(all(baseline.Pass), ...
    'Frozen local inputs failed to reproduce canonical lesion EDSS rows.');

membership_table = struct2table(membership_rows);
membership_table = sortrows(membership_table, ...
    {'ThresholdLabel','Metric','Model','SubjectID'});
assert(~any(startsWith(membership_table.SubjectID,'sub-C')), ...
    'A control entered the lesion EDSS sensitivity cohort.');

writetable(result_table,fullfile(output_dir, ...
    'lesion_edss_threshold_model_results.csv'));
writetable(baseline,fullfile(qa_dir, ...
    'lesion_edss_baseline_reproduction.csv'));
writetable(membership_table,fullfile(qa_dir, ...
    'lesion_edss_subject_membership.csv'));

report = struct();
report.schema_version = 1;
report.status = 'PASS';
report.scope = ['Eight unadjusted lesion-load EDSS models, canonical >3 ' ...
    'versus sensitivity >=4'];
report.canonical_rows_reproduced = height(baseline);
report.canonical_rows_all_pass = all(baseline.Pass);
report.canonical_rule = 'EDSS > 3.0';
report.sensitivity_rule = 'EDSS >= 4.0';
report.N = unique(result_table.N);
report.canonical_positive_N = unique(result_table.PositiveN( ...
    result_table.ThresholdLabel == "canonical_gt3"));
report.sensitivity_positive_N = unique(result_table.PositiveN( ...
    result_table.ThresholdLabel == "sensitivity_ge4"));
report.model_count_per_threshold = 8;
report.bootstrap = ['Legacy unstratified 2000-resample bootstrap mean ' ...
    'and percentile interval; rng(42) reset per model.'];
report.cohort_note = ['RR047 is retained only to hold the current canonical ' ...
    'Supplementary Table 5 cohort fixed while isolating the threshold change.'];
report.read_only_inputs = struct('clinical',clinical_file, ...
    'tract_lesion',tract_file,'classical_lesion',classical_file, ...
    'canonical_models',canonical_file);
write_json(fullfile(qa_dir,'lesion_edss_ge4_report.json'),report);
fprintf('Lesion EDSS >=4 isolated sensitivity PASS.\n');
end


function value = is_text_scalar(value)
value = ischar(value) || (isstring(value) && isscalar(value));
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


function res = run_model_binary(X,y,n_boot,alpha)
rng(42);
X_z = zscore(X);
lambda_grid = logspace(-6,6,50);
cv = fitclinear(X_z,y,'Learner','logistic', ...
    'Regularization','ridge','Lambda',lambda_grid,'KFold',10);
[~,best] = min(kfoldLoss(cv));
mdl_r = fitclinear(X_z,y,'Learner','logistic', ...
    'Regularization','ridge','Lambda',lambda_grid(best));
LP = X_z*mdl_r.Beta + mdl_r.Bias;
warning_state = warning;
cleanup = onCleanup(@() warning(warning_state));
warning('off','all');
mdl = fitglm(LP,y,'Distribution','binomial');
score = LP;
coef = mdl.Coefficients.Estimate(2);
se = mdl.Coefficients.SE(2);
p = 2*normcdf(-abs(coef/se));
AIC = mdl.ModelCriterion.AIC;
[auc_mean,auc_ci] = auc_ci_bootstrap(y,score,n_boot,alpha);
res = struct('p',p,'AIC',AIC,'AUC_mean',auc_mean, ...
    'AUC_CI',auc_ci,'lambda',lambda_grid(best));
end


function [mean_auc,ci] = auc_ci_bootstrap(y,score,n_boot,alpha)
y = y(:);
score = score(:);
rng(42,'twister');
aucs = bootstrp(n_boot,@(yy,ss) perf_auc(yy,ss),y,score);
mean_auc = mean(aucs);
ci = quantile(aucs,[alpha/2,1-alpha/2]);
end


function value = perf_auc(y,score)
[~,~,~,value] = perfcurve(y,score,1);
end


function [X,y,ids,edss] = get_complete_cases(T,preds,measure,threshold,use_geq)
X_raw = table2array(T(:,preds));
y_raw = double(T.(measure));
valid = all(isfinite(X_raw),2) & isfinite(y_raw);
X = X_raw(valid,:);
ids = string(T.SubjectID(valid));
edss = y_raw(valid);
if use_geq, y = double(edss >= threshold); else, y = double(edss > threshold); end
assert(numel(unique(y)) == 2,'Binary outcome must contain both classes.');
end


function [Xn,Xc,y,ids,edss] = get_joint_cases( ...
        Tn,Tc,nrf_preds,crf_preds,measure,threshold,use_geq)
[~,ia,ib] = intersect(Tn.SubjectID,Tc.SubjectID,'stable');
Tn = Tn(ia,:);
Tc = Tc(ib,:);
Xn_raw = table2array(Tn(:,nrf_preds));
Xc_raw = table2array(Tc(:,crf_preds));
y_raw = double(Tn.(measure));
valid = all(isfinite(Xn_raw),2) & all(isfinite(Xc_raw),2) & ...
    isfinite(y_raw);
Xn = Xn_raw(valid,:);
Xc = Xc_raw(valid,:);
ids = string(Tn.SubjectID(valid));
edss = y_raw(valid);
if use_geq, y = double(edss >= threshold); else, y = double(edss > threshold); end
assert(numel(unique(y)) == 2,'Binary outcome must contain both classes.');
end


function preds = build_nrf_preds(groups,metric)
preds = {};
for gi = 1:numel(groups)
    preds{end+1} = [groups{gi} 'L' metric]; %#ok<AGROW>
    preds{end+1} = [groups{gi} 'R' metric]; %#ok<AGROW>
end
end


function preds = build_crf_preds(regions,metric)
preds = {};
for ri = 1:numel(regions)
    preds{end+1} = [regions{ri} metric]; %#ok<AGROW>
end
end


function validate_predictors(T,preds,label)
missing = preds(~ismember(preds,T.Properties.VariableNames));
assert(isempty(missing),'Missing %s predictors: %s',label,strjoin(missing,', '));
end


function rows = append_result(rows,threshold_label,positive_rule, ...
        metric,model,ids,y,res)
r = struct();
r.ThresholdLabel = threshold_label;
r.PositiveRule = positive_rule;
r.Outcome = 'EDSS';
r.Metric = metric;
r.Model = model;
r.Adjusted = 0;
r.Demographics = 'none';
r.N = numel(y);
r.PositiveN = sum(y);
r.PerformanceName = 'AUC';
r.Performance = res.AUC_mean;
r.CI_lower = res.AUC_CI(1);
r.CI_upper = res.AUC_CI(2);
r.p_raw = res.p;
r.AIC = res.AIC;
r.Lambda = res.lambda;
r.SubjectIDSetSHA256 = sha256_text(strjoin(sort(ids),newline));
rows = append_struct(rows,r);
end


function rows = append_membership(rows,threshold_label,metric,model,ids,edss,y)
for i = 1:numel(ids)
    r = struct();
    r.ThresholdLabel = threshold_label;
    r.Metric = metric;
    r.Model = model;
    r.SubjectID = char(ids(i));
    r.EDSS = edss(i);
    r.Positive = y(i);
    rows = append_struct(rows,r);
end
end


function T = sort_result_table(T)
threshold_order = double(T.ThresholdLabel == "sensitivity_ge4") + 1;
metric_names = ["LN","LV","Lnorm"];
metric_order = zeros(height(T),1);
for i = 1:numel(metric_names)
    metric_order(T.Metric == metric_names(i)) = i;
end
model_order = zeros(height(T),1);
model_order(T.Model == "Whole Brain") = 1;
model_order(T.Model == "Tract-based") = 2;
model_order(T.Model == "Classical-region") = 3;
T.ThresholdOrder = threshold_order;
T.MetricOrder = metric_order;
T.ModelOrder = model_order;
T = sortrows(T,{'ThresholdOrder','MetricOrder','ModelOrder'});
T.ThresholdOrder = [];
T.MetricOrder = [];
T.ModelOrder = [];
end


function baseline = validate_baseline(results,canonical_file)
canonical = readtable(canonical_file,'TextType','string');
canonical = canonical(canonical.Outcome == "EDSS",:);
observed = results(results.ThresholdLabel == "canonical_gt3",:);
assert(height(canonical) == 8 && height(observed) == 8, ...
    'Expected eight canonical and eight reproduced EDSS lesion rows.');
checks = struct([]);
for i = 1:height(canonical)
    key = observed.Metric == canonical.Metric(i) & ...
        observed.Model == canonical.Model(i);
    assert(nnz(key) == 1,'Canonical lesion key is not unique.');
    row = observed(key,:);
    c = struct();
    c.Metric = char(canonical.Metric(i));
    c.Model = char(canonical.Model(i));
    c.NMatch = row.N == canonical.N(i);
    c.PositiveNMatch = row.PositiveN == canonical.PositiveN(i);
    c.PerformanceMatch = abs(row.Performance-canonical.Performance(i)) < 5e-4;
    c.CILowerMatch = abs(row.CI_lower-canonical.CI_lower(i)) < 5e-4;
    c.CIUpperMatch = abs(row.CI_upper-canonical.CI_upper(i)) < 5e-4;
    c.RawPMatch = abs(row.p_raw-canonical.p_raw(i)) < 1e-12;
    c.AICMatch = abs(row.AIC-canonical.AIC(i)) < 5e-3;
    c.ReproducedPerformance = row.Performance;
    c.CanonicalPerformance = canonical.Performance(i);
    c.ReproducedP = row.p_raw;
    c.CanonicalP = canonical.p_raw(i);
    c.ReproducedAIC = row.AIC;
    c.CanonicalAIC = canonical.AIC(i);
    c.Pass = c.NMatch && c.PositiveNMatch && c.PerformanceMatch && ...
        c.CILowerMatch && c.CIUpperMatch && c.RawPMatch && c.AICMatch;
    checks = append_struct(checks,c);
end
baseline = struct2table(checks);
end


function values = append_struct(values,item)
if isempty(values), values = item; else, values(end+1) = item; end
end


function value = ternary(condition,a,b)
if condition, value = a; else, value = b; end
end


function h = sha256_text(value)
md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
bytes = unicode2native(char(value),'UTF-8');
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
