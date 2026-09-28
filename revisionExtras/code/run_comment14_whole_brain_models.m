function report = run_comment14_whole_brain_models(varargin)
%RUN_COMMENT14_WHOLE_BRAIN_MODELS Fit Comment 14 whole-brain qMRI models.
%
% This reviewer-only analysis uses one whole-brain parenchymal mean for each
% of T1, MTR, FA, and MD. Models contain the imaging metric alone, use no
% demographic covariates, exclude controls, and reuse the joint
% tract/classical complete-case cohorts used by the manuscript comparisons.

p = inputParser;
addParameter(p, 'ClinicalFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'WholeBrainFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'MetricsDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ClassicalFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'QADir', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

code_dir = fileparts(mfilename('fullpath'));
package_root = fileparts(code_dir);
clinical_file = char(p.Results.ClinicalFile);
whole_brain_file = char(p.Results.WholeBrainFile);
metrics_dir = char(p.Results.MetricsDir);
classical_file = char(p.Results.ClassicalFile);
output_dir = char(p.Results.OutputDir);
qa_dir = char(p.Results.QADir);
if isempty(classical_file)
    classical_file = fullfile(metrics_dir, 'ClassicalRegionAllMetrics.csv');
end
if isempty(output_dir)
    output_dir = fullfile(package_root, 'work', 'comment14_manual', 'model_results');
end
if isempty(qa_dir)
    qa_dir = fullfile(package_root, 'work', 'comment14_manual', 'qa');
end
output_dir = require_inside_package(output_dir, package_root, 'OutputDir');
qa_dir = require_inside_package(qa_dir, package_root, 'QADir');
assert(isfile(clinical_file), 'Clinical input missing: %s', clinical_file);
assert(isfile(whole_brain_file), 'Whole-brain input missing: %s', whole_brain_file);
assert(isfolder(metrics_dir), 'Metric input directory missing: %s', metrics_dir);
assert(isfile(classical_file), 'Classical-region input missing: %s', classical_file);
if ~isfolder(output_dir), mkdir(output_dir); end
if ~isfolder(qa_dir), mkdir(qa_dir); end

rng(42, 'twister');
n_boot = 2000;
alpha_ci = 0.05;
metrics = {'T1','MTR','FA','MD'};
bin_measures = {'EDSS','MSPro'};
bin_labels = {'EDSS','MSPro'};
bin_thresholds = [4,1];
bin_use_geq = [true,false];
cont_measures = {'T25FW','x9HPTD','x9HPTND'};
cont_labels = {'T25FW','9HPT-D','9HPT-ND'};
cont_transforms = {'log','log','log'};
expected_binary_n = [80,80,79,79; 79,79,78,78];
expected_binary_positive = [31,31,31,31; 39,39,39,39];
expected_continuous_n = [76,76,75,75; 79,79,78,78; 79,79,78,78];
tract_preds = build_tract_predictors();

T_clin = readtable(clinical_file);
assert(height(T_clin) == 132 && ...
    numel(unique(string(T_clin.SubjectID))) == 132, ...
    'Clinical table must contain 132 unique frozen subjects.');

T_wb = readtable(whole_brain_file);
assert(height(T_wb) == 132 && numel(unique(string(T_wb.SubjectID))) == 132, ...
    'Whole-brain table must contain 132 unique frozen imaging subjects.');
assert(~any(string(T_wb.SubjectID) == "sub-RR047"), ...
    'RR047 must remain excluded from whole-brain inputs.');
assert(all(isfinite(T_wb.T1)) && all(isfinite(T_wb.FA)) && ...
    all(isfinite(T_wb.MD)), 'T1/FA/MD must be complete for the frozen 132.');
expected_mtr_missing = ismember(string(T_wb.SubjectID), ["sub-C035","sub-RR013"]);
assert(sum(expected_mtr_missing) == 2 && ...
    all(~isfinite(T_wb.MTR(expected_mtr_missing))) && ...
    all(isfinite(T_wb.MTR(~expected_mtr_missing))), ...
    'MTR must be missing for C035/RR013 only.');

T_base = outerjoin(T_clin, T_wb, 'Keys', 'SubjectID', 'MergeKeys', true);
T_classical = readtable(classical_file);
classical_ids = string(T_classical.SubjectID);
assert(height(T_classical) == 87 && numel(unique(classical_ids)) == 87, ...
    'Classical-region table must contain 87 unique MS participants.');
assert(~any(startsWith(classical_ids, 'sub-C')), ...
    'Classical-region cohort must not contain controls.');
model_rows = struct([]);
tract_hashes = struct();

for mi = 1:numel(metrics)
    metric = metrics{mi};
    tract_file = fullfile(metrics_dir, sprintf('GroupTract%s_All.xlsx', metric));
    assert(isfile(tract_file), 'Tract input missing: %s', tract_file);
    T_tract = readtable(tract_file);
    tract_ids = string(T_tract.SubjectID);
    assert(height(T_tract) == 132 && numel(unique(tract_ids)) == 132, ...
        'Tract table for %s must contain 132 unique subjects.', metric);
    assert(isempty(setdiff(tract_preds, T_tract.Properties.VariableNames)), ...
        'Tract table for %s is missing required predictors.', metric);
    classical_preds = build_classical_predictors(metric);
    assert(isempty(setdiff(classical_preds, T_classical.Properties.VariableNames)), ...
        'Classical-region table is missing required %s predictors.', metric);
    tract_hashes.(metric) = sha256_file(tract_file);

    for oi = 1:numel(bin_measures)
        outcome = bin_labels{oi};
        measure = bin_measures{oi};
        [X,y,ids] = get_matched_cases(T_base, T_tract, T_classical, ...
            tract_preds, classical_preds, metric, measure, true, ...
            bin_thresholds(oi), bin_use_geq(oi), 'raw');
        assert(numel(y) == expected_binary_n(oi,mi) && ...
            numel(unique(y)) == 2 && ...
            sum(y) == expected_binary_positive(oi,mi), ...
            'Binary cohort mismatch for %s/%s.', outcome, metric);
        rng(42, 'twister');
        [LP,lambda] = fit_ridge_binary(zscore(X), y);
        result = evaluate_binary(LP, y, [], n_boot, alpha_ci);
        model_rows = append_struct(model_rows, make_model_row( ...
            outcome, measure, metric, false, {}, result, lambda, ...
            true, y, ids));
    end
    for oi = 1:numel(cont_measures)
        outcome = cont_labels{oi};
        measure = cont_measures{oi};
        transform = cont_transforms{oi};
        [X,y,ids] = get_matched_cases(T_base, T_tract, T_classical, ...
            tract_preds, classical_preds, metric, measure, false, ...
            NaN, false, transform);
        assert(numel(y) == expected_continuous_n(oi,mi), ...
            'Continuous cohort mismatch for %s/%s.', outcome, metric);
        rng(42, 'twister');
        result = run_model_continuous(X, y, []);
        model_rows = append_struct(model_rows, make_model_row( ...
            outcome, measure, metric, false, {}, result, ...
            result.lambda, false, [], ids));
    end
end

model_table = order_model_rows(struct2table(model_rows));
assert(height(model_table) == 20, 'Expected 20 whole-brain model rows.');
assert(nnz(model_table.Adjusted) == 0, 'All models must be unadjusted.');
assert(nnz(strcmp(model_table.Status, 'COMPLETED')) == 20, ...
    'All whole-brain models must complete.');
output_file = fullfile(output_dir, 'whole_brain_qmri_models.csv');
writetable(model_table, output_file);

report = struct();
report.schema_version = 2;
report.status = 'PASS';
report.analysis_role = 'reviewer-only Comment 14 whole-brain qMRI models';
report.edss_rule = 'EDSS >= 4.0 (revisionExtras sensitivity definition)';
report.metrics = metrics;
report.outcomes = [bin_labels cont_labels];
report.model_specification = 'one whole-brain imaging metric; no covariates';
report.cohort_policy = ['joint tract+classical complete-case MS cohorts: ' ...
    'Supplementary Table 3 cohorts for EDSS/MSPro and analogous matched ' ...
    'manuscript cohorts for T25FW/9HPT'];
report.controls_included = false;
report.covariates = 'none';
report.rng_seed = 42;
report.binary_lambda_grid = 'logspace(-6,6,50)';
report.continuous_lambda_grid = 'logspace(-6,6,100)';
report.cv_folds = 10;
report.bootstrap_repetitions = n_boot;
report.multiplicity = 'raw p only; outside the manuscript multiplicity family';
report.whole_brain_input_sha256 = sha256_file(whole_brain_file);
report.clinical_input_sha256 = sha256_file(clinical_file);
report.classical_input_sha256 = sha256_file(classical_file);
report.tract_input_sha256 = tract_hashes;
report.output_rows = height(model_table);
report.adjusted_rows = 0;
report.raw_p_lt_0_05 = nnz(model_table.p_raw < 0.05);
write_json(fullfile(qa_dir, 'comment14_whole_brain_report.json'), report);
fprintf('Comment 14 whole-brain modeling complete: %d rows.\n', height(model_table));
end


function [X,y,ids] = get_matched_cases(T_base,T_tract,T_classical, ...
        tract_preds,classical_preds,metric,measure,is_binary,threshold, ...
        use_geq,transform)
[ids_bt,ib,it] = intersect(string(T_base.SubjectID), ...
    string(T_tract.SubjectID),'stable');
B = T_base(ib,:);
R = T_tract(it,:);
[ids_all,ij,ic] = intersect(ids_bt,string(T_classical.SubjectID),'stable');
B = B(ij,:);
R = R(ij,:);
C = T_classical(ic,:);
assert(isequal(string(B.SubjectID),string(R.SubjectID)) && ...
    isequal(string(B.SubjectID),string(C.SubjectID)), ...
    'Matched cohort alignment failed for %s.', metric);

Xraw = double(B.(metric));
yraw = double(B.(measure));
Xtract = table2array(R(:,tract_preds));
Xclassical = table2array(C(:,classical_preds));
valid = ~startsWith(ids_all,'sub-C') & isfinite(Xraw) & Xraw~=0 & ...
    isfinite(yraw) & all(isfinite(Xtract) & Xtract~=0,2) & ...
    all(isfinite(Xclassical) & Xclassical~=0,2);
if ~is_binary && strcmp(transform,'log')
    valid = valid & yraw>0;
end
X = Xraw(valid,:);
ids = ids_all(valid);
if is_binary
    if use_geq
        y = double(yraw(valid)>=threshold);
    else
        y = double(yraw(valid)>threshold);
    end
elseif strcmp(transform,'log')
    y = log(yraw(valid));
else
    y = yraw(valid);
end
end


function preds = build_tract_predictors()
groups = {'Association','Cerebellar','Occipitoparietal','ProjectionBrainstem'};
preds = {};
for i = 1:numel(groups)
    group = groups{i};
    preds = [preds, {[group 'L_Tail'],[group 'R_Tail'], ...
        [group 'L_NAWM'],[group 'R_NAWM']}]; %#ok<AGROW>
end
end


function preds = build_classical_predictors(metric)
regions = {'periventricular','juxtacortical','infratentorial','deepwhitematter'};
preds = cellfun(@(region) [region metric], regions, 'UniformOutput', false);
end


function [LP,lambda] = fit_ridge_binary(Xz,y)
grid = logspace(-6,6,50);
cv = fitclinear(Xz,y,'Learner','logistic','Regularization','ridge', ...
    'Lambda',grid,'KFold',10);
[~,index] = min(kfoldLoss(cv)); lambda = grid(index);
model = fitclinear(Xz,y,'Learner','logistic','Regularization','ridge', ...
    'Lambda',lambda);
LP = Xz*model.Beta+model.Bias;
end


function result = evaluate_binary(LP,y,Demo,n_boot,alpha_ci)
warning('off','all');
if isempty(Demo)
    model = fitglm(LP,y,'Distribution','binomial'); score = LP;
else
    model = fitglm([LP,Demo],y,'Distribution','binomial');
    score = model.Fitted.LinearPredictor;
end
warning('on','all');
coef = model.Coefficients.Estimate(2); se = model.Coefficients.SE(2);
p = 2*normcdf(-abs(coef/se));
rng(42,'twister');
aucs = bootstrp(n_boot,@(yy,ss) perf_auc(yy,ss),y(:),score(:));
result = struct('n',numel(y),'p',p,'AIC',model.ModelCriterion.AIC, ...
    'Performance',mean(aucs),'CI',[quantile(aucs,alpha_ci/2), ...
    quantile(aucs,1-alpha_ci/2)]);
end


function value = perf_auc(y,score)
[~,~,~,value] = perfcurve(y,score,1);
end


function result = run_model_continuous(X,y,Demo)
Xz = zscore(X); grid = logspace(-6,6,100);
cv = fitrlinear(Xz,y,'Learner','leastsquares','Regularization','ridge', ...
    'Lambda',grid,'KFold',10);
[~,index] = min(kfoldLoss(cv)); lambda = grid(index);
ridge = fitrlinear(Xz,y,'Learner','leastsquares','Regularization','ridge', ...
    'Lambda',lambda);
LP = Xz*ridge.Beta+ridge.Bias;
if isempty(Demo)
    model = fitlm(LP,y); R2 = model.Rsquared.Ordinary;
else
    model = fitlm([LP,Demo],y); R2 = fitlm(model.Fitted,y).Rsquared.Ordinary;
end
coef = model.Coefficients.Estimate(2); se = model.Coefficients.SE(2);
t = coef/se; p = 2*tcdf(-abs(t),model.DFE);
if p==0, p = 2*normcdf(-abs(t)); end
result = struct('n',numel(y),'p',p,'AIC',model.ModelCriterion.AIC, ...
    'Performance',R2,'CI',[NaN,NaN],'lambda',lambda);
end


function row = make_model_row(outcome,source,metric,adjusted,demos, ...
        result,lambda,is_binary,y,ids)
row = struct('AnalysisScope','matched_manuscript_case','Outcome',outcome, ...
    'SourceField',source,'Metric',metric,'Model','Whole brain', ...
    'Adjusted',adjusted,'Demographics','','CasePolicy', ...
    'joint tract+classical complete-case MS cohort','N',result.n, ...
    'PositiveN',NaN,'NegativeN',NaN,'PerformanceName','', ...
    'Performance',result.Performance,'CI_lower',result.CI(1), ...
    'CI_upper',result.CI(2),'Lambda',lambda,'p_raw',result.p, ...
    'AIC',result.AIC,'SubjectIDSetSHA256',hash_ids(ids), ...
    'Status','COMPLETED','Blocker','');
row.Demographics = strjoin(demos,'+');
if isempty(row.Demographics), row.Demographics = 'none'; end
if is_binary
    row.PositiveN = sum(y); row.NegativeN = numel(y)-sum(y);
    row.PerformanceName = 'apparent AUC';
else
    row.PerformanceName = 'apparent R2';
end
end


function T = order_model_rows(T)
metric_order = categorical(string(T.Metric),{'T1','MTR','FA','MD'},'Ordinal',true);
outcome_order = categorical(string(T.Outcome), ...
    {'EDSS','MSPro','T25FW','9HPT-D','9HPT-ND'}, ...
    'Ordinal',true);
T.MetricOrderTemp = metric_order; T.OutcomeOrderTemp = outcome_order;
T = sortrows(T,{'MetricOrderTemp','OutcomeOrderTemp','Adjusted'});
T(:,{'MetricOrderTemp','OutcomeOrderTemp'}) = [];
end


function values = append_struct(values,item)
if isempty(values), values=item; else, values(end+1)=item; end
end


function value = hash_ids(ids)
ids = sort(string(ids(:)));
bytes = unicode2native(char(strjoin(ids,newline)),'UTF-8');
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(bytes); raw = typecast(engine.digest(),'uint8');
value = lower(reshape(dec2hex(raw,2).',1,[]));
end


function digest = sha256_file(path)
fid = fopen(path,'rb'); assert(fid>=0,'Could not hash %s.',path);
cleanup = onCleanup(@()fclose(fid)); bytes = fread(fid,Inf,'*uint8');
engine = java.security.MessageDigest.getInstance('SHA-256'); engine.update(bytes);
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end


function write_json(path,value)
fid = fopen(path,'w'); assert(fid>=0,'Could not write %s.',path);
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end


function value = require_inside_package(path,package_root,label)
value = char(java.io.File(path).getCanonicalPath());
root = char(java.io.File(package_root).getCanonicalPath());
assert(startsWith([value filesep],[root filesep]), ...
    '%s must remain inside revisionExtras: %s',label,value);
end
