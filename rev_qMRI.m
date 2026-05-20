function results = rev_qMRI()
% REV_QMRI  Binary + continuous ridge regression for qMRI metrics (T1, MTR, FA, MD).
%
% Two model families per metric:
%   Tract-based — 16 predictors: 4 groups × L/R × {Tail, NAWM}
%   Classic     —  4 predictors: 4 icometrix regions × metric
%
% BINARY outcomes  : EDSS (> 3), MSPro (> 1)   — logistic ridge, 50-point lambda grid
% CONTINUOUS outcomes: T25FW, 9HPT-D, 9HPT-ND  — linear ridge, 100-point lambda grid,
%                      log-transformed outcome
%
% Controls (sub-C) excluded for binary; retained for continuous.
%   Binary   : joint complete cases — both Tract-based and Classic share the same N.
%              Matches paper_multi_bi_ridge.m (check_predictors='both').
%   Continuous: independent complete cases — each model has its own N.
%              Matches paper_multi_continuous_ridge.m (check_predictors='model1/2').
%
% No-demo vs +Demo rows:
%   Each is an INDEPENDENT ridge regression run, matching the original approach of
%   calling paper_multi_bi_ridge_all('demographics',{}) for no-demo and
%   paper_multi_bi_ridge_all() for +Demo. N may differ between the two rows if
%   any subjects have missing Age/Gender data.
%
% Zeros excluded from qMRI predictor columns (intensity metrics).
% rng(42) reset once per NRF+Classic pair (before the NRF call) — matches the
% original paper_multi_bi_ridge.m / paper_multi_continuous_ridge.m behavior where
% rng(42) is set at function entry and NRF runs before Classic in sequence.
% Demographics: Age + Gender (no DurationOfDisease —
%   matches paper_multi_continuous_ridge_all.m default).
%
% Usage:
%   rev_qMRI

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

qmri_metrics  = {'T1', 'MTR', 'FA', 'MD'};

% NRF predictor names — column names in GroupTract{metric}_All.xlsx.
% The metric is encoded in which file we load, NOT in the column name.
% Groups verified from actual file headers.
tract_groups_full = {'Association','Cerebellar','Occipitoparietal','ProjectionBrainstem'};

% CRF predictor names — columns in icometrixLesionStats.xlsx.
% Pattern: {region}{metric}, e.g. periventricularT1.
% Verified from actual file headers.
ico_regions = {'periventricular','juxtacortical','infratentorial','deepwhitematter'};

% NRF predictor list (same for every metric — metric encoded in the file)
nrf_preds_base = build_nrf_preds(tract_groups_full);   % 16 names

demo_cols_bin  = {'Age','GenderNum','DurationOfDisease'};  % matches paper_multi_bi_ridge.m default
demo_cols_cont = {'Age','GenderNum'};                       % matches paper_multi_continuous_ridge_all.m default
n_boot    = 2000;
alpha_ci  = 0.05;

%% ── Load clinical data ────────────────────────────────────────────────────
fprintf('Loading clinical data...\n');
T_clin = readtable(fullfile(base_dir, 'clinicalScore.xlsx'));

% Gender encoding: M=1, F=0, empty string→NaN.
% NaN matches load_ridge_data.m lines 227-228: cell columns with empty values are
% treated as missing and excluded.  Subjects with blank Gender would otherwise
% receive GenderNum=0 (Female) and survive the demo filter incorrectly.
if iscell(T_clin.Gender) || isstring(T_clin.Gender)
    gender_empty = cellfun(@(x) isempty(x), cellstr(T_clin.Gender));
    T_clin.GenderNum = double(strcmp(T_clin.Gender, 'M'));
    T_clin.GenderNum(gender_empty) = NaN;
else
    T_clin.GenderNum = double(T_clin.Gender == 1);
end

%% ── Load CRF data (one shared file for all metrics) ──────────────────────
% icometrixLesionStats.xlsx contains columns for all qMRI metrics:
%   periventricularT1, periventricularMTR, ... deepwhitematterMD, etc.
fprintf('Loading CRF data (icometrixLesionStats.xlsx)...\n');
T_ico      = readtable(fullfile(base_dir, 'icometrixLesionStats.xlsx'));
T_crf_full = outerjoin(T_clin, T_ico, 'Keys','SubjectID','MergeKeys',true);
T_crf_ms   = T_crf_full(~startsWith(T_crf_full.SubjectID, 'sub-C'), :);

%% ── Load NRF data per metric (GroupTract{metric}_All.xlsx) ───────────────
% Each metric has its own file; column names are metric-independent.
fprintf('Loading NRF data per metric...\n');
T_nrf_full = struct();   % controls retained  (for continuous)
T_nrf_ms   = struct();   % controls removed   (for binary)
for mi = 1:numel(qmri_metrics)
    metric   = qmri_metrics{mi};
    nrf_file = fullfile(base_dir, sprintf('GroupTract%s_All.xlsx', metric));
    if ~exist(nrf_file, 'file')
        error('NRF file not found: %s', nrf_file);
    end
    T_tract = readtable(nrf_file);
    T_full  = outerjoin(T_clin, T_tract, 'Keys','SubjectID','MergeKeys',true);
    T_nrf_full.(metric) = T_full;
    T_nrf_ms.(metric)   = T_full(~startsWith(T_full.SubjectID, 'sub-C'), :);
    fprintf('  %s: %d subjects total, %d MS patients\n', metric, ...
        height(T_full), height(T_nrf_ms.(metric)));
end

% Validate NRF predictor names against the first metric's table
nrf_preds_base = validate_preds(nrf_preds_base, T_nrf_full.(qmri_metrics{1}), 'NRF-qMRI');

%% ── Output ────────────────────────────────────────────────────────────────
results  = struct();
all_rows = {};   % each model → 2 rows: no-demo and +Demo (independent runs)

%% ════════════════════════════════════════════════════════════════════════
%  BINARY OUTCOMES
%  Joint complete cases: subjects need finite, non-zero data in BOTH NRF
%  and CRF predictor columns (matches paper_multi_bi_ridge.m logic).
%  No-demo and +Demo are independent ridge runs (different N if any subjects
%  have missing Age/Gender).
%% ════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n  BINARY OUTCOMES\n%s\n', repmat('=',1,65), repmat('=',1,65));

for ci = 1:numel(bin_measures)
    measure   = bin_measures{ci};
    threshold = bin_thresh(ci);
    fprintf('\n>>> %s (threshold > %d)\n', measure, threshold);

    for mi = 1:numel(qmri_metrics)
        metric    = qmri_metrics{mi};
        nrf_preds = validate_preds(nrf_preds_base, T_nrf_ms.(metric), sprintf('NRF-%s', metric));
        crf_preds = validate_preds(build_crf_preds(ico_regions, metric), ...
                                   T_crf_ms, sprintf('CRF-%s', metric));
        if isempty(nrf_preds) || isempty(crf_preds), continue; end

        name_nrf = sprintf('Tract-based %s', metric);
        name_crf = sprintf('Classic %s', metric);

        % ── No-demo: no demographic filtering (demographics={} in original)
        [Xn_nd, Xc_nd, y_nd, ~, n_nd] = get_joint_cases( ...
            T_nrf_ms.(metric), T_crf_ms, nrf_preds, crf_preds, ...
            measure, {}, true, threshold, use_geq, true);

        lambda_grid_bin = logspace(-6, 6, 50);   % 50-point grid, matches paper_multi_bi_ridge.m

        % ── No-demo pair ──────────────────────────────────────────────────
        % rng(42) once → NRF fitclinear → Classic fitclinear (inherits state).
        % eval_binary (fitglm + AUC) is called AFTER both ridge fits so that
        % auc_ci_bootstrap's internal rng(42,'twister') reset does not corrupt
        % the state needed by Classic fitclinear.
        % Matches paper_multi_bi_ridge.m: rng at top, NRF fitclinear, Classic
        % fitclinear, then AUC for all four models.
        if n_nd >= 20
            rng(42);
            LP_nrf_nd = fit_ridge_lp_bin(zscore(Xn_nd), y_nd, lambda_grid_bin, ...
                                         sprintf('%s|%s|nd', measure, name_nrf));
            LP_crf_nd = fit_ridge_lp_bin(zscore(Xc_nd), y_nd, lambda_grid_bin, ...
                                         sprintf('%s|%s|nd', measure, name_crf));

            res = eval_binary(LP_nrf_nd, y_nd, [], n_nd, n_boot, alpha_ci, ...
                              sprintf('%s|%s|nodemo', measure, name_nrf));
            results.(sprintf('%s_NRF_%s_nd', measure, metric)) = res;
            all_rows = append_bin_row(all_rows, measure, name_nrf, res);

            res = eval_binary(LP_crf_nd, y_nd, [], n_nd, n_boot, alpha_ci, ...
                              sprintf('%s|%s|nodemo', measure, name_crf));
            results.(sprintf('%s_CRF_%s_nd', measure, metric)) = res;
            all_rows = append_bin_row(all_rows, measure, name_crf, res);
        else
            warn_skip(measure, sprintf('%s/%s no-demo', name_nrf, name_crf), n_nd);
        end

        % ── With-demo pair ────────────────────────────────────────────────
        [Xn_wd, Xc_wd, y_wd, Demo_wd, n_wd] = get_joint_cases( ...
            T_nrf_ms.(metric), T_crf_ms, nrf_preds, crf_preds, ...
            measure, demo_cols_bin, true, threshold, use_geq, true);

        if n_wd >= 20
            rng(42);
            LP_nrf_wd = fit_ridge_lp_bin(zscore(Xn_wd), y_wd, lambda_grid_bin, ...
                                         sprintf('%s|%s|wd', measure, name_nrf));
            LP_crf_wd = fit_ridge_lp_bin(zscore(Xc_wd), y_wd, lambda_grid_bin, ...
                                         sprintf('%s|%s|wd', measure, name_crf));

            res = eval_binary(LP_nrf_wd, y_wd, Demo_wd, n_wd, n_boot, alpha_ci, ...
                              sprintf('%s|%s|demo', measure, name_nrf));
            results.(sprintf('%s_NRF_%s_wd', measure, metric)) = res;
            all_rows = append_bin_row(all_rows, measure, [name_nrf ' (+Demo)'], res);

            res = eval_binary(LP_crf_wd, y_wd, Demo_wd, n_wd, n_boot, alpha_ci, ...
                              sprintf('%s|%s|demo', measure, name_crf));
            results.(sprintf('%s_CRF_%s_wd', measure, metric)) = res;
            all_rows = append_bin_row(all_rows, measure, [name_crf ' (+Demo)'], res);
        else
            warn_skip(measure, sprintf('%s/%s +demo', name_nrf, name_crf), n_wd);
        end
    end
end

%% ════════════════════════════════════════════════════════════════════════
%  CONTINUOUS OUTCOMES
%  Independent complete cases: NRF and CRF filtered separately → different N.
%  Matches paper_multi_continuous_ridge.m (check_predictors='model1/model2').
%  No-demo and +Demo are also independent ridge runs.
%% ════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n  CONTINUOUS OUTCOMES (log-transformed)\n%s\n', repmat('=',1,65), repmat('=',1,65));

for ci = 1:numel(cont_measures)
    measure = cont_measures{ci};
    clabel  = cont_labels{ci};
    fprintf('\n>>> %s\n', clabel);

    for mi = 1:numel(qmri_metrics)
        metric    = qmri_metrics{mi};
        nrf_preds = validate_preds(nrf_preds_base, T_nrf_full.(metric), sprintf('NRF-%s', metric));
        crf_preds = validate_preds(build_crf_preds(ico_regions, metric), ...
                                   T_crf_full, sprintf('CRF-%s', metric));

        name_nrf = sprintf('Tract-based %s', metric);
        name_crf = sprintf('Classic %s', metric);

        % Load all four case sets upfront (data loading is deterministic — no rng)
        [Xn_nd, yn_nd, ~,    nn_nd] = get_complete_cases( ...
            T_nrf_full.(metric), nrf_preds, measure, {}, false, [], [], true);
        [Xn_wd, yn_wd, Dn_wd, nn_wd] = get_complete_cases( ...
            T_nrf_full.(metric), nrf_preds, measure, demo_cols_cont, false, [], [], true);
        [Xc_nd, yc_nd, ~,    nc_nd] = get_complete_cases( ...
            T_crf_full, crf_preds, measure, {}, false, [], [], true);
        [Xc_wd, yc_wd, Dc_wd, nc_wd] = get_complete_cases( ...
            T_crf_full, crf_preds, measure, demo_cols_cont, false, [], [], true);

        % No-demo pair: rng(42) → NRF fitrlinear → Classic fitrlinear.
        % Matches paper_multi_continuous_ridge.m (rng(42) at top, NRF then Classic).
        rng(42);
        if ~isempty(nrf_preds) && nn_nd >= 20
            res = run_model_cont(Xn_nd, yn_nd, [], nn_nd, sprintf('%s|%s|nd', clabel, name_nrf));
            results.(sprintf('%s_NRF_%s_nd', measure, metric)) = res;
            all_rows = append_cont_row(all_rows, clabel, name_nrf, res);
        elseif ~isempty(nrf_preds)
            warn_skip(clabel, [name_nrf ' no-demo'], nn_nd);
        end
        if ~isempty(crf_preds) && nc_nd >= 20
            res = run_model_cont(Xc_nd, yc_nd, [], nc_nd, sprintf('%s|%s|nd', clabel, name_crf));
            results.(sprintf('%s_CRF_%s_nd', measure, metric)) = res;
            all_rows = append_cont_row(all_rows, clabel, name_crf, res);
        elseif ~isempty(crf_preds)
            warn_skip(clabel, [name_crf ' no-demo'], nc_nd);
        end

        % +Demo pair: fresh rng(42) → NRF fitrlinear → Classic fitrlinear.
        rng(42);
        if ~isempty(nrf_preds) && nn_wd >= 20
            res = run_model_cont(Xn_wd, yn_wd, Dn_wd, nn_wd, sprintf('%s|%s|wd', clabel, name_nrf));
            results.(sprintf('%s_NRF_%s_wd', measure, metric)) = res;
            all_rows = append_cont_row(all_rows, clabel, [name_nrf ' (+Demo)'], res);
        elseif ~isempty(nrf_preds)
            warn_skip(clabel, [name_nrf ' +demo'], nn_wd);
        end
        if ~isempty(crf_preds) && nc_wd >= 20
            res = run_model_cont(Xc_wd, yc_wd, Dc_wd, nc_wd, sprintf('%s|%s|wd', clabel, name_crf));
            results.(sprintf('%s_CRF_%s_wd', measure, metric)) = res;
            all_rows = append_cont_row(all_rows, clabel, [name_crf ' (+Demo)'], res);
        elseif ~isempty(crf_preds)
            warn_skip(clabel, [name_crf ' +demo'], nc_wd);
        end
    end
end

%% ── Print & save ──────────────────────────────────────────────────────────
print_combined_table(all_rows);
out_file = fullfile(script_dir, 'rev_qMRI_results.xlsx');
save_excel(out_file, all_rows);
fprintf('\nResults saved to: %s\n', out_file);

end  % ── main ──────────────────────────────────────────────────────────────


%% ══════════════════════════════════════════════════════════════════════════
%  ANALYSIS FUNCTIONS
%% ══════════════════════════════════════════════════════════════════════════

function LP = fit_ridge_lp_bin(X_z, y, lambda_grid, tag)
% Phase 1: fit ridge logistic and return the linear predictor LP.
%
% rng is NOT reset here.  The caller sets rng(42) once before the NRF call so
% that NRF consumes rng state first; the Classic call then inherits that state —
% matching paper_multi_bi_ridge.m (rng at top, NRF fitclinear, then Classic fitclinear).
% auc_ci_bootstrap is intentionally NOT called here; it resets rng(42,'twister')
% internally and would corrupt the state needed by Classic fitclinear.

    cv    = fitclinear(X_z, y, 'Learner','logistic', ...
                'Regularization','ridge', 'Lambda',lambda_grid, 'KFold',10);
    [~, bi] = min(kfoldLoss(cv));
    mdl_r   = fitclinear(X_z, y, 'Learner','logistic', ...
                'Regularization','ridge', 'Lambda',lambda_grid(bi));
    LP = X_z * mdl_r.Beta + mdl_r.Bias;
    fprintf('  [ridge] %s  lambda=%.2e\n', tag, lambda_grid(bi));
end


function res = eval_binary(LP, y, Demo, n, n_boot, alpha_ci, tag)
% Phase 2: fit GLM and compute AUC.  rng state does not matter here because
% auc_ci_bootstrap resets rng(42,'twister') internally every call.

    has_demo = ~isempty(Demo);
    n_pos    = sum(y);
    fprintf('\n  [binary%s] %s  N=%d  pos=%d (%.1f%%)\n', ...
        ternary(has_demo,'+demo',''), tag, n, n_pos, 100*n_pos/n);

    warning('off','all');
    if has_demo
        mdl   = fitglm([LP, Demo], y, 'Distribution','binomial');
        score = mdl.Fitted.LinearPredictor;
    else
        mdl   = fitglm(LP, y, 'Distribution','binomial');
        score = LP;
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
% y is already log-transformed.
%
% rng is NOT reset here — caller resets rng(42) once before the NRF call so
% that the Classic call inherits the rng state after NRF, matching paper_multi_continuous_ridge.m.

    has_demo = ~isempty(Demo);
    fprintf('\n  [cont%s] %s  N=%d  mean=%.3f SD=%.3f\n', ...
        ternary(has_demo,'+demo',''), tag, n, mean(y), std(y));

    X_z         = zscore(X);
    lambda_grid = logspace(-6, 6, 100);  % 100-point grid — matches paper_multi_continuous_ridge.m

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
% Columns: Outcome | Model | N (pos) | p | AUC [95%CI] | AIC | _type
    rows{end+1, 1} = outcome;
    rows{end,   2} = model_name;
    rows{end,   3} = sprintf('%d (%d)', res.n, res.n_pos);
    rows{end,   4} = res.p;   % exact numeric — needed for BH correction
    rows{end,   5} = sprintf('%.3f [%.3f, %.3f]', res.AUC_mean, res.AUC_CI(1), res.AUC_CI(2));
    rows{end,   6} = sprintf('%.2f', res.AIC);
    rows{end,   7} = 'bin';
end

function rows = append_cont_row(rows, outcome, model_name, res)
% Columns: Outcome | Model | N | p | R² | AIC | _type
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
% Return filtered X, y, Demo for a single table.
% demo_cols={} → no demographic filtering (matches demographics={} original run).
% excl_pred_zeros: exclude rows where any predictor == 0 (intensity metrics).

    X_raw = table2array(T(:, pred_names));
    y_raw = T.(measure);

    nd = numel(demo_cols);
    Dmat = zeros(height(T), nd);
    for di = 1:nd, Dmat(:,di) = T.(demo_cols{di}); end

    pred_ok = all(isfinite(X_raw), 2);
    if excl_pred_zeros
        pred_ok = pred_ok & all(X_raw ~= 0, 2);
    end
    % all(isfinite(Dmat), 2) returns true for all rows when Dmat is Nx0
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
% Joint complete cases: subjects need finite (+ non-zero) data in BOTH model sets.
% demo_cols={} → no demographic filtering (matches demographics={} original run).
% Matches paper_multi_bi_ridge.m check_predictors='both' with same metric.

    [~, ia, ib] = intersect(T_nrf.SubjectID, T_crf.SubjectID, 'stable');
    T_n = T_nrf(ia, :);
    T_c = T_crf(ib, :);

    X_nrf_raw = table2array(T_n(:, nrf_preds));
    X_crf_raw = table2array(T_c(:, crf_preds));
    y_raw     = T_n.(measure);

    nd = numel(demo_cols);
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
    Demo  = Dmat(valid, :);   % Nx0 when demo_cols={}
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

function preds = build_nrf_preds(tract_groups_full)
% 16 predictors: 4 groups × L/R × {Tail, NAWM}.
% Column names verified from GroupTract{metric}_All.xlsx headers.
    preds = {};
    for gi = 1:numel(tract_groups_full)
        g = tract_groups_full{gi};
        preds{end+1} = [g 'L_Tail'];  %#ok<AGROW>
        preds{end+1} = [g 'R_Tail'];  %#ok<AGROW>
        preds{end+1} = [g 'L_NAWM']; %#ok<AGROW>
        preds{end+1} = [g 'R_NAWM']; %#ok<AGROW>
    end
end

function preds = build_crf_preds(ico_regions, metric)
% 4 predictors: {region}{metric}, e.g. periventricularT1.
% Column names verified from icometrixLesionStats.xlsx headers.
    preds = {};
    for ri = 1:numel(ico_regions)
        preds{end+1} = [ico_regions{ri} metric]; %#ok<AGROW>
    end
end

function preds = validate_preds(preds, T, tag)
    missing = preds(~ismember(preds, T.Properties.VariableNames));
    if ~isempty(missing)
        warning('rev_qMRI:%s — missing columns: %s', tag, strjoin(missing, ', '));
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
    rng(42, 'twister');   % matches paper_multi_bi_ridge.m auc_ci_bootstrap
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
    warning('rev_qMRI: only %d cases for %s — %s. Skipping.', n, outcome, model);
end

function s = pval_str(p)
    if p < 0.001, s = '<0.001'; else, s = sprintf('%.4f', p); end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
