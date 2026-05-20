function results = rev_lesionload()
% REV_LESIONLOAD  Binary + continuous ridge regression for lesion load metrics.
%
% Three model families:
%   WB          — whole-brain (WBLN, WBLV)           
%   Tract-based — tract-based (LN, LV, Lnorm)       
%   Classic     — icometrix regional (LN, LV, Lnorm) 
%
% BINARY outcomes  : EDSS (> 3), MSPro (> 1)   — logistic ridge, 50-point lambda grid
% CONTINUOUS outcomes: T25FW, 9HPT-D, 9HPT-ND  — linear ridge, 100-point lambda grid,
%                      log-transformed outcome
%
% Controls (sub-C) excluded for binary; retained for continuous.
%   Binary   : joint complete cases — Tract-based and Classic share same N.
%              Matches paper_multi_bi_ridge.m (check_predictors='both').
%   Continuous: joint complete cases for Tract-based and Classic (LN/LV/Lnorm).
%
% No-demo vs +Demo rows:
%   Each is an INDEPENDENT ridge regression run, matching the original approach of
%   calling paper_multi_bi_ridge_all('demographics',{}) for no-demo and
%   paper_multi_bi_ridge_all() for +Demo. N may differ between the two rows if
%   any subjects have missing Age/Gender data.
%
% Demographics: Age + Gender (no DurationOfDisease —
%   matches paper_multi_continuous_ridge_all.m default).
% rng(42) reset at start of each model run.
%
% Output: single combined table (binary then continuous) printed and saved
%         to rev_lesionload_results.xlsx in the same folder as this script.
%
% Usage:
%   rev_lesionload

rng(42, 'twister');

%% ── Paths ─────────────────────────────────────────────────────────────────
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, 'helpers')));
cfg      = load_config();
base_dir = fullfile(cfg.bidsDir, cfg.statsDir);

%% ── Config ────────────────────────────────────────────────────────────────
bin_measures = {'EDSS',  'MSPro'};
bin_thresh   = [3,        1     ];
use_geq      = false;   % y = score > threshold  (matches paper_multi_bi_ridge.m)

cont_measures = {'T25FW',  'x9HPTD',  'x9HPTND'};
cont_labels   = {'T25FW',  '9HPT-D',  '9HPT-ND'};

tract_groups   = {'Association','Cerebellar','Occipitoparietal','PB'};
ico_regions    = {'periventricular','juxtacortical','infratentorial','deepwhitematter'};
lesion_metrics = {'LN','LV','Lnorm'};

demo_cols_bin  = {'Age','GenderNum','DurationOfDisease'};  % matches paper_multi_bi_ridge.m default
demo_cols_cont = {'Age','GenderNum'};                       % matches paper_multi_continuous_ridge_all.m default
n_boot = 2000;   alpha_ci = 0.05;

%% ── Load data ─────────────────────────────────────────────────────────────
fprintf('Loading data...\n');
T_clin  = readtable(fullfile(base_dir, 'clinicalScore.xlsx'));
T_tract = readtable(fullfile(base_dir, 'GroupTractLesionLoad.xlsx'));
T_ico   = readtable(fullfile(base_dir, 'icometrixLesionLoad.xlsx'));

% Gender encoding: M=1, F=0, empty string→NaN.
% NaN matches load_ridge_data.m lines 227-228: empty-string cell values are
% treated as missing — subjects with blank Gender are excluded from demo filter.
if iscell(T_clin.Gender) || isstring(T_clin.Gender)
    gender_empty = cellfun(@(x) isempty(x), cellstr(T_clin.Gender));
    T_clin.GenderNum = double(strcmp(T_clin.Gender, 'M'));
    T_clin.GenderNum(gender_empty) = NaN;
else
    T_clin.GenderNum = double(T_clin.Gender == 1);
end

T_nrf = outerjoin(T_clin, T_tract, 'Keys','SubjectID','MergeKeys',true);
T_crf = outerjoin(T_clin, T_ico,   'Keys','SubjectID','MergeKeys',true);

ctrl_nrf = startsWith(T_nrf.SubjectID, 'sub-C');
ctrl_crf = startsWith(T_crf.SubjectID, 'sub-C');
T_nrf_ms = T_nrf(~ctrl_nrf, :);
T_crf_ms = T_crf(~ctrl_crf, :);
fprintf('Controls removed: %d (tract table), %d (ico table). MS patients: %d\n', ...
    sum(ctrl_nrf), sum(ctrl_crf), height(T_nrf_ms));

%% ── Output ────────────────────────────────────────────────────────────────
results  = struct();
all_rows = {};   % each model → 2 rows: no-demo and +Demo (independent runs)

%% ════════════════════════════════════════════════════════════════════════
%  BINARY OUTCOMES
%% ════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n  BINARY OUTCOMES\n%s\n', repmat('=',1,65), repmat('=',1,65));

for ci = 1:numel(bin_measures)
    measure   = bin_measures{ci};
    threshold = bin_thresh(ci);
    fprintf('\n>>> %s (threshold > %d)\n', measure, threshold);

    % ── WB models (single predictor, no joint filtering needed)
    for wi = 1:2
        wb_col   = ternary(wi==1,'WBLN','WBLV');
        name_wb  = wb_col;

        [X_nd, y_nd, ~,     n_nd] = get_complete_cases(T_nrf_ms, {wb_col}, measure, {}, true, threshold, use_geq, false);
        [X_wd, y_wd, D_wd, n_wd] = get_complete_cases(T_nrf_ms, {wb_col}, measure, demo_cols_bin, true, threshold, use_geq, false);

        if n_nd >= 20
            res = run_model_binary(X_nd, y_nd, [], n_nd, n_boot, alpha_ci, sprintf('%s|%s|nd', measure, name_wb));
            results.(sprintf('%s_%s_nd', measure, wb_col)) = res;
            all_rows = append_bin_row(all_rows, measure, name_wb, res);
        else, warn_skip(measure, [name_wb ' no-demo'], n_nd); end

        if n_wd >= 20
            res = run_model_binary(X_wd, y_wd, D_wd, n_wd, n_boot, alpha_ci, sprintf('%s|%s|wd', measure, name_wb));
            results.(sprintf('%s_%s_wd', measure, wb_col)) = res;
            all_rows = append_bin_row(all_rows, measure, [name_wb ' (+Demo)'], res);
        else, warn_skip(measure, [name_wb ' +demo'], n_wd); end
    end

    % ── Tract-based and Classic (joint complete cases)
    for mi = 1:numel(lesion_metrics)
        metric    = lesion_metrics{mi};
        nrf_preds = validate_preds(build_nrf_preds(tract_groups, metric), T_nrf_ms, sprintf('Tract-%s', metric));
        crf_preds = validate_preds(build_crf_preds(ico_regions,  metric), T_crf_ms, sprintf('Classic-%s', metric));
        if isempty(nrf_preds) || isempty(crf_preds), continue; end

        name_nrf = sprintf('Tract-based %s', metric);
        name_crf = sprintf('Classic %s', metric);

        % No-demo: no demographic filtering
        [Xn_nd, Xc_nd, y_nd, ~,     n_nd] = get_joint_cases(T_nrf_ms, T_crf_ms, nrf_preds, crf_preds, measure, {}, true, threshold, use_geq, false);
        % With-demo: filter on Age + Gender + DurationOfDisease (matches paper_multi_bi_ridge.m)
        [Xn_wd, Xc_wd, y_wd, D_wd, n_wd] = get_joint_cases(T_nrf_ms, T_crf_ms, nrf_preds, crf_preds, measure, demo_cols_bin, true, threshold, use_geq, false);

        if n_nd >= 20
            res = run_model_binary(Xn_nd, y_nd, [], n_nd, n_boot, alpha_ci, sprintf('%s|%s|nd', measure, name_nrf));
            results.(sprintf('%s_NRF_%s_nd', measure, metric)) = res;
            all_rows = append_bin_row(all_rows, measure, name_nrf, res);

            res = run_model_binary(Xc_nd, y_nd, [], n_nd, n_boot, alpha_ci, sprintf('%s|%s|nd', measure, name_crf));
            results.(sprintf('%s_CRF_%s_nd', measure, metric)) = res;
            all_rows = append_bin_row(all_rows, measure, name_crf, res);
        else, warn_skip(measure, sprintf('%s/%s no-demo', name_nrf, name_crf), n_nd); end

        if n_wd >= 20
            res = run_model_binary(Xn_wd, y_wd, D_wd, n_wd, n_boot, alpha_ci, sprintf('%s|%s|wd', measure, name_nrf));
            results.(sprintf('%s_NRF_%s_wd', measure, metric)) = res;
            all_rows = append_bin_row(all_rows, measure, [name_nrf ' (+Demo)'], res);

            res = run_model_binary(Xc_wd, y_wd, D_wd, n_wd, n_boot, alpha_ci, sprintf('%s|%s|wd', measure, name_crf));
            results.(sprintf('%s_CRF_%s_wd', measure, metric)) = res;
            all_rows = append_bin_row(all_rows, measure, [name_crf ' (+Demo)'], res);
        else, warn_skip(measure, sprintf('%s/%s +demo', name_nrf, name_crf), n_wd); end
    end
end

%% ════════════════════════════════════════════════════════════════════════
%  CONTINUOUS OUTCOMES
%% ════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n  CONTINUOUS OUTCOMES (log-transformed)\n%s\n', repmat('=',1,65), repmat('=',1,65));

for ci = 1:numel(cont_measures)
    measure = cont_measures{ci};
    clabel  = cont_labels{ci};
    fprintf('\n>>> %s\n', clabel);

    % ── WB models (controls retained for continuous)
    for wi = 1:2
        wb_col  = ternary(wi==1,'WBLN','WBLV');
        name_wb = wb_col;

        [X_nd, y_nd, ~,     n_nd] = get_complete_cases(T_nrf, {wb_col}, measure, {}, false, [], [], false);
        [X_wd, y_wd, D_wd, n_wd] = get_complete_cases(T_nrf, {wb_col}, measure, demo_cols_cont, false, [], [], false);

        if n_nd >= 20
            res = run_model_cont(X_nd, y_nd, [], n_nd, sprintf('%s|%s|nd', clabel, name_wb));
            results.(sprintf('%s_%s_nd', measure, wb_col)) = res;
            all_rows = append_cont_row(all_rows, clabel, name_wb, res);
        else, warn_skip(clabel, [name_wb ' no-demo'], n_nd); end

        if n_wd >= 20
            res = run_model_cont(X_wd, y_wd, D_wd, n_wd, sprintf('%s|%s|wd', clabel, name_wb));
            results.(sprintf('%s_%s_wd', measure, wb_col)) = res;
            all_rows = append_cont_row(all_rows, clabel, [name_wb ' (+Demo)'], res);
        else, warn_skip(clabel, [name_wb ' +demo'], n_wd); end
    end

    % ── Tract-based and Classic (joint complete cases, controls retained)
    for mi = 1:numel(lesion_metrics)
        metric    = lesion_metrics{mi};
        nrf_preds = validate_preds(build_nrf_preds(tract_groups, metric), T_nrf, sprintf('Tract-%s', metric));
        crf_preds = validate_preds(build_crf_preds(ico_regions,  metric), T_crf, sprintf('Classic-%s', metric));
        if isempty(nrf_preds) || isempty(crf_preds), continue; end

        name_nrf = sprintf('Tract-based %s', metric);
        name_crf = sprintf('Classic %s', metric);

        [Xn_nd, Xc_nd, y_nd, ~,     n_nd] = get_joint_cases(T_nrf, T_crf, nrf_preds, crf_preds, measure, {}, false, [], [], false);
        [Xn_wd, Xc_wd, y_wd, D_wd, n_wd] = get_joint_cases(T_nrf, T_crf, nrf_preds, crf_preds, measure, demo_cols_cont, false, [], [], false);

        if n_nd >= 20
            res = run_model_cont(Xn_nd, y_nd, [], n_nd, sprintf('%s|%s|nd', clabel, name_nrf));
            results.(sprintf('%s_NRF_%s_nd', measure, metric)) = res;
            all_rows = append_cont_row(all_rows, clabel, name_nrf, res);

            res = run_model_cont(Xc_nd, y_nd, [], n_nd, sprintf('%s|%s|nd', clabel, name_crf));
            results.(sprintf('%s_CRF_%s_nd', measure, metric)) = res;
            all_rows = append_cont_row(all_rows, clabel, name_crf, res);
        else, warn_skip(clabel, sprintf('%s/%s no-demo', name_nrf, name_crf), n_nd); end

        if n_wd >= 20
            res = run_model_cont(Xn_wd, y_wd, D_wd, n_wd, sprintf('%s|%s|wd', clabel, name_nrf));
            results.(sprintf('%s_NRF_%s_wd', measure, metric)) = res;
            all_rows = append_cont_row(all_rows, clabel, [name_nrf ' (+Demo)'], res);

            res = run_model_cont(Xc_wd, y_wd, D_wd, n_wd, sprintf('%s|%s|wd', clabel, name_crf));
            results.(sprintf('%s_CRF_%s_wd', measure, metric)) = res;
            all_rows = append_cont_row(all_rows, clabel, [name_crf ' (+Demo)'], res);
        else, warn_skip(clabel, sprintf('%s/%s +demo', name_nrf, name_crf), n_wd); end
    end
end

%% ── Print & save ──────────────────────────────────────────────────────────
print_combined_table(all_rows);
out_file = fullfile(script_dir, 'rev_lesionload_results.xlsx');
save_excel(out_file, all_rows);
fprintf('\nResults saved to: %s\n', out_file);

end  % ── main ──────────────────────────────────────────────────────────────


%% ══════════════════════════════════════════════════════════════════════════
%  ANALYSIS FUNCTIONS
%% ══════════════════════════════════════════════════════════════════════════

function res = run_model_binary(X, y, Demo, n, n_boot, alpha_ci, tag)
% Ridge logistic, one run.  Demo=[] → no-demo model; Demo provided → +Demo model.
% rng(42) reset per call — matches paper_multi_bi_ridge.m (fresh call per metric).

    rng(42);
    has_demo = ~isempty(Demo);
    n_pos    = sum(y);
    fprintf('\n  [binary%s] %s  N=%d  pos=%d (%.1f%%)\n', ...
        ternary(has_demo,'+demo',''), tag, n, n_pos, 100*n_pos/n);

    X_z         = zscore(X);
    lambda_grid = logspace(-6, 6, 50);

    cv    = fitclinear(X_z, y, 'Learner','logistic', ...
                'Regularization','ridge', 'Lambda',lambda_grid, 'KFold',10);
    [~, bi] = min(kfoldLoss(cv));
    mdl_r   = fitclinear(X_z, y, 'Learner','logistic', ...
                'Regularization','ridge', 'Lambda',lambda_grid(bi));
    LP1 = X_z * mdl_r.Beta + mdl_r.Bias;
    fprintf('  lambda=%.2e\n', lambda_grid(bi));

    warning('off','all');
    if has_demo
        mdl   = fitglm([LP1, Demo], y, 'Distribution','binomial');
        score = mdl.Fitted.LinearPredictor;
    else
        mdl   = fitglm(LP1, y, 'Distribution','binomial');
        score = LP1;
    end
    warning('on','all');

    [~, ~, p, AIC, ~] = extract_binary_stats(mdl);
    [auc_m, auc_ci]   = auc_ci_bootstrap(y, score, n_boot, alpha_ci);
    fprintf('  p=%s  AUC=%.3f [%.3f,%.3f]  AIC=%.2f\n', ...
        pval_str(p), auc_m, auc_ci(1), auc_ci(2), AIC);

    res = struct('n',n, 'n_pos',n_pos, 'p',p, 'AIC',AIC, 'AUC_mean',auc_m, 'AUC_CI',auc_ci);
end


function res = run_model_cont(X, y, Demo, n, tag)
% Ridge linear, one run.  Demo=[] → no-demo model; Demo provided → +Demo model.
% y is already log-transformed.  rng(42) reset per call.

    rng(42);
    has_demo = ~isempty(Demo);
    fprintf('\n  [cont%s] %s  N=%d  mean=%.3f SD=%.3f\n', ...
        ternary(has_demo,'+demo',''), tag, n, mean(y), std(y));

    X_z         = zscore(X);
    lambda_grid = logspace(-6, 6, 100);

    cv    = fitrlinear(X_z, y, 'Learner','leastsquares', ...
                'Regularization','ridge', 'Lambda',lambda_grid, 'KFold',10);
    [~, bi] = min(kfoldLoss(cv));
    mdl_r   = fitrlinear(X_z, y, 'Learner','leastsquares', ...
                'Regularization','ridge', 'Lambda',lambda_grid(bi));
    LP1 = X_z * mdl_r.Beta + mdl_r.Bias;
    fprintf('  lambda=%.2e\n', lambda_grid(bi));

    if has_demo
        mdl  = fitlm([LP1, Demo], y);
        LP11 = mdl.Fitted;
        [~, ~, p, AIC, ~, ~] = extract_linear_stats(mdl);
        R2   = fitlm(LP11, y).Rsquared.Ordinary;
    else
        mdl  = fitlm(LP1, y);
        [~, ~, p, AIC, R2, ~] = extract_linear_stats(mdl);
    end
    fprintf('  p=%s  R²=%.3f  AIC=%.2f\n', pval_str(p), R2, AIC);

    res = struct('n',n, 'p',p, 'AIC',AIC, 'R2',R2);
end


%% ══════════════════════════════════════════════════════════════════════════
%  ROW BUILDERS  (one row per call)
%% ══════════════════════════════════════════════════════════════════════════

function rows = append_bin_row(rows, outcome, model_name, res)
    rows{end+1, 1} = outcome;
    rows{end,   2} = model_name;
    rows{end,   3} = sprintf('%d (%d)', res.n, res.n_pos);
    rows{end,   4} = res.p;   % exact numeric — needed for BH correction
    rows{end,   5} = sprintf('%.3f [%.3f, %.3f]', res.AUC_mean, res.AUC_CI(1), res.AUC_CI(2));
    rows{end,   6} = sprintf('%.2f', res.AIC);
    rows{end,   7} = 'bin';
end

function rows = append_cont_row(rows, outcome, model_name, res)
    rows{end+1, 1} = outcome;
    rows{end,   2} = model_name;
    rows{end,   3} = sprintf('%d', res.n);
    rows{end,   4} = res.p;   % exact numeric — needed for BH correction
    rows{end,   5} = sprintf('%.3f', res.R2);
    rows{end,   6} = sprintf('%.2f', res.AIC);
    rows{end,   7} = 'cont';
end


%% ══════════════════════════════════════════════════════════════════════════
%  DATA LOADING
%% ══════════════════════════════════════════════════════════════════════════

function [X, y, Demo, n] = get_complete_cases(T, pred_names, measure, ...
                                               demo_cols, is_binary, threshold, use_geq, ...
                                               excl_pred_zeros)
% demo_cols={} → no demographic filtering (matches demographics={} original run).
% excl_pred_zeros: exclude rows where any predictor == 0.
%   false for lesion metrics (zeros = no lesions, valid).
%   true for intensity metrics (zeros indicate missing data).

    X_raw = table2array(T(:, pred_names));
    y_raw = T.(measure);

    nd   = numel(demo_cols);
    Dmat = zeros(height(T), nd);
    for di = 1:nd, Dmat(:,di) = T.(demo_cols{di}); end

    pred_ok = all(isfinite(X_raw), 2);
    if excl_pred_zeros
        pred_ok = pred_ok & all(X_raw ~= 0, 2);
    end
    % all(isfinite(Nx0 matrix), 2) = true(N,1) — no filtering when demo_cols={}
    demo_ok = all(isfinite(Dmat), 2);

    if is_binary
        valid = pred_ok & isfinite(y_raw) & demo_ok;
    else
        valid = pred_ok & isfinite(y_raw) & y_raw > 0 & demo_ok;
    end

    X    = X_raw(valid, :);
    Demo = Dmat(valid, :);   % Nx0 when demo_cols={}; isempty(Demo) → true
    n    = sum(valid);
    y_v  = y_raw(valid);
    if is_binary
        if use_geq, y = double(y_v >= threshold); else, y = double(y_v > threshold); end
    else
        y = log(y_v);
    end
end


function [X_nrf, X_crf, y, Demo, n] = get_joint_cases(T_nrf, T_crf, ...
        nrf_preds, crf_preds, measure, demo_cols, is_binary, threshold, use_geq, ...
        excl_pred_zeros)
% Joint complete cases: subjects need finite data in BOTH NRF and CRF predictors.
% demo_cols={} → no demographic filtering.

    [~, ia, ib] = intersect(T_nrf.SubjectID, T_crf.SubjectID, 'stable');
    T_n = T_nrf(ia, :);
    T_c = T_crf(ib, :);

    X_nrf_raw = table2array(T_n(:, nrf_preds));
    X_crf_raw = table2array(T_c(:, crf_preds));
    y_raw     = T_n.(measure);

    nd   = numel(demo_cols);
    Dmat = zeros(height(T_n), nd);
    for di = 1:nd, Dmat(:,di) = T_n.(demo_cols{di}); end

    pred_ok = all(isfinite(X_nrf_raw), 2) & all(isfinite(X_crf_raw), 2);
    if excl_pred_zeros
        pred_ok = pred_ok & all(X_nrf_raw ~= 0, 2) & all(X_crf_raw ~= 0, 2);
    end
    demo_ok = all(isfinite(Dmat), 2);

    if is_binary
        valid = pred_ok & isfinite(y_raw) & demo_ok;
    else
        valid = pred_ok & isfinite(y_raw) & y_raw > 0 & demo_ok;
    end

    X_nrf = X_nrf_raw(valid, :);
    X_crf = X_crf_raw(valid, :);
    Demo  = Dmat(valid, :);
    n     = sum(valid);
    y_v   = y_raw(valid);
    if is_binary
        if use_geq, y = double(y_v >= threshold); else, y = double(y_v > threshold); end
    else
        y = log(y_v);
    end
end


%% ══════════════════════════════════════════════════════════════════════════
%  PREDICTOR BUILDERS
%% ══════════════════════════════════════════════════════════════════════════

function preds = build_nrf_preds(tract_groups, metric)
% 8 predictors: 4 groups × L/R.  Column format: {Group}L{metric}, e.g. AssociationLLN.
    preds = {};
    for gi = 1:numel(tract_groups)
        preds{end+1} = [tract_groups{gi} 'L' metric]; %#ok<AGROW>
        preds{end+1} = [tract_groups{gi} 'R' metric]; %#ok<AGROW>
    end
end

function preds = build_crf_preds(ico_regions, metric)
% 4 predictors: {region}{metric}, e.g. periventricularLN.
    preds = {};
    for ri = 1:numel(ico_regions)
        preds{end+1} = [ico_regions{ri} metric]; %#ok<AGROW>
    end
end

function preds = validate_preds(preds, T, tag)
    missing = preds(~ismember(preds, T.Properties.VariableNames));
    if ~isempty(missing)
        warning('rev_lesionload:%s — missing: %s', tag, strjoin(missing,', '));
        preds = preds(ismember(preds, T.Properties.VariableNames));
    end
end


%% ══════════════════════════════════════════════════════════════════════════
%  STAT EXTRACTORS
%% ══════════════════════════════════════════════════════════════════════════

function [OR, CI, p, AIC, cal] = extract_binary_stats(mdl)
% Manual Wald z-test for p — avoids MATLAB fitglm NaN-p bug (binomial+Dispersion=1).
% Uses lower-tail normcdf to avoid upper-tail underflow (1-normcdf(z)=0 for z>~8.2).
    coef = mdl.Coefficients.Estimate(2);
    se   = mdl.Coefficients.SE(2);
    OR   = exp(coef);
    CI   = exp(coef + [-1 1]*1.96*se);
    p    = 2 * normcdf(-abs(coef / se));   % stable lower-tail form
    AIC  = mdl.ModelCriterion.AIC;
    cal  = coef;
end

function [beta, CI, p, AIC, R2, cal] = extract_linear_stats(mdl)
    coef = mdl.Coefficients.Estimate(2);
    se   = mdl.Coefficients.SE(2);
    df   = mdl.DFE;
    t    = coef / se;
    beta = coef;
    CI   = coef + [-1 1]*1.96*se;
    % tcdf(-|t|, df) is the stable form but still underflows for very large |t|
    % (tcdf uses betainc internally). Fall back to normcdf (accurate for df>=30).
    p = 2 * tcdf(-abs(t), df);
    if p == 0
        p = 2 * normcdf(-abs(t));
    end
    AIC  = mdl.ModelCriterion.AIC;
    R2   = mdl.Rsquared.Ordinary;
    cal  = coef;
end


%% ══════════════════════════════════════════════════════════════════════════
%  AUC BOOTSTRAP
%% ══════════════════════════════════════════════════════════════════════════

function [AUC_mean, AUC_CI] = auc_ci_bootstrap(y, scores, n_boot, alpha_ci)
    y = y(:);  scores = scores(:);
    rng(42, 'twister');
    AUCs     = bootstrp(n_boot, @(yy,ss) perf_auc(yy,ss), y, scores);
    AUC_mean = mean(AUCs);
    AUC_CI   = quantile(AUCs, [alpha_ci/2, 1-alpha_ci/2]);
end

function a = perf_auc(y, s)
    [~,~,~,a] = perfcurve(y, s, 1);
end


%% ══════════════════════════════════════════════════════════════════════════
%  PRINT COMBINED TABLE
%% ══════════════════════════════════════════════════════════════════════════

function print_combined_table(rows)
    if isempty(rows), return; end

    w   = [10, 26, 12, 8, 28, 10];
    hdr = {'Outcome', 'Model', 'N', 'p', 'AUC [95%CI] / R²', 'AIC'};
    fmt = '';
    for i = 1:numel(w), fmt = [fmt sprintf('%%-%ds  ', w(i))]; end %#ok<AGROW>
    fmt = [fmt '\n'];
    sep = repmat('-', 1, sum(w) + 2*numel(w));

    fprintf('\n\n%s\n  COMBINED RESULTS TABLE\n%s\n', repmat('=',1,65), repmat('=',1,65));
    fprintf(fmt, hdr{:});
    fprintf('%s\n', sep);

    prev_type    = '';
    prev_outcome = '';
    for r = 1:size(rows, 1)
        cur_type    = rows{r, 7};
        cur_outcome = rows{r, 1};

        if ~isempty(prev_type) && ~strcmp(cur_type, prev_type)
            fprintf('\n');
            fprintf(fmt, hdr{:});
            fprintf('%s\n', sep);
        end

        if strcmp(cur_type, prev_type) && ~isempty(prev_outcome) && ...
                ~strcmp(cur_outcome, prev_outcome)
            fprintf('\n');
        end

        fprintf(fmt, rows{r, 1:6});
        prev_type    = cur_type;
        prev_outcome = cur_outcome;
    end
    fprintf('\n');
end


%% ══════════════════════════════════════════════════════════════════════════
%  SAVE (single Excel sheet)
%% ══════════════════════════════════════════════════════════════════════════

function save_excel(out_file, rows)
    if isempty(rows), return; end
    col_names = {'Outcome','Model','N','p','AUC_CI_or_R2','AIC'};

    % Split no-demo rows (sheet "Results") from +Demo rows (sheet "Demo")
    is_demo  = cellfun(@(m) contains(m, '(+Demo)'), rows(:, 2));
    rows_nd  = rows(~is_demo, :);
    rows_d   = rows( is_demo, :);

    if ~isempty(rows_nd)
        writetable(cell2table(rows_nd(:, 1:6), 'VariableNames', col_names), ...
                   out_file, 'Sheet', 'Results');
    end
    if ~isempty(rows_d)
        writetable(cell2table(rows_d(:, 1:6), 'VariableNames', col_names), ...
                   out_file, 'Sheet', 'Demo');
    end
end


%% ══════════════════════════════════════════════════════════════════════════
%  UTILITIES
%% ══════════════════════════════════════════════════════════════════════════

function warn_skip(outcome, model, n)
    warning('rev_lesionload: only %d cases for %s — %s. Skipping.', n, outcome, model);
end

function s = pval_str(p)
    if p < 0.001, s = '<0.001'; else, s = sprintf('%.4f', p); end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
