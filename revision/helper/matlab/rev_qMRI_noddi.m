function report = rev_qMRI_noddi(varargin)
% REV_QMRI_NODDI  Source-derived qMRI analysis, with a unified all-metric mode.
%
% This file was copied from the public MS100 rev_qMRI.m and minimally adapted:
%   1. input/output directories are name-value arguments;
%   2. NDI, ODI, and ISOVF (file label FWF) are supported alongside the
%      original T1, MTR, FA, and MD metrics;
%   3. the same 16-predictor tract and four-predictor classical branches run;
%   4. raw SDMT and signed raw MSFC-SDMT are appended to the five original
%      outcomes; and
%   5. results are saved in the reviewer-package canonical CSV schema.
%
% The ridge fitting, lambda grids, RNG seed, full-sample refit, apparent
% AUC/R2, demographic covariates, and case rules below retain the mechanics
% of rev_qMRI.m. Exact-zero NODDI estimates are valid and are not discarded.
%
% BINARY outcomes  : EDSS (> 3), MSPro (> 1)   — logistic ridge, 50-point lambda grid
% CONTINUOUS outcomes: log(T25FW), log(9HPT-D), log(9HPT-ND), raw SDMT,
%                      and signed raw MSFC-SDMT — linear ridge, 100-point grid
%
% Controls (sub-C) are excluded for binary and retained for continuous.
% Binary models use joint tract/classical complete cases. Continuous
% models use independent available cases; exact matched-subject refits are
% produced separately for Figure 6/S2. FA/MD are refit first with the original
% legacy case rules as an equivalence gate.
%
% No-demo vs +Demo rows:
%   Each is an INDEPENDENT ridge regression run, matching the original approach of
%   calling paper_multi_bi_ridge_all('demographics',{}) for no-demo and
%   paper_multi_bi_ridge_all() for +Demo. N may differ between the two rows if
%   any subjects have missing Age/Gender data.
%
% Demographics: Age + Gender for continuous, and Age + Gender + disease
% duration for binary, exactly as in rev_qMRI.m.
%
% Unified publication analysis (preferred):
%   run_all_qmri_manuscript_models(...)
%
% Internal NODDI-only compatibility mode (must be requested explicitly):
%   rev_qMRI_noddi('statsDir', statsDir, 'metricsDir', metricsDir, ...
%                  'outputDir', outputDir, ...
%                  'classicalFile', classicalFile, ...
%                  'analysisScope', 'noddi-only')

p = inputParser;
addParameter(p, 'statsDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'metricsDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'classicalFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'outputDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'nBoot', 2000, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'globalFamilySize', 98, ...
    @(x) isnumeric(x) && isscalar(x) && x == 98);
addParameter(p, 'analysisScope', 'all-metrics', ...
    @(x) any(strcmpi(char(x), {'noddi-only','all-metrics'})));
addParameter(p, 'referenceValidationMode', 'frozen', ...
    @(x) any(strcmpi(char(x), {'frozen','candidate-migration'})));
parse(p, varargin{:});

stats_dir = char(p.Results.statsDir);
metrics_dir = char(p.Results.metricsDir);
classical_file = char(p.Results.classicalFile);
output_dir = char(p.Results.outputDir);
analysis_scope = lower(char(p.Results.analysisScope));
reference_validation_mode = lower(char(p.Results.referenceValidationMode));
all_metrics_mode = strcmp(analysis_scope, 'all-metrics');
assert(~isempty(stats_dir), 'statsDir is required.');
assert(~isempty(metrics_dir), 'metricsDir is required.');
assert(~isempty(output_dir), 'outputDir is required.');
if isempty(classical_file)
    classical_file = fullfile(metrics_dir, 'ClassicalRegionNODDI_All.csv');
end
if ~isfolder(output_dir), mkdir(output_dir); end

rng(42, 'twister');

%% ── Paths ─────────────────────────────────────────────────────────────────
%% ── Config ────────────────────────────────────────────────────────────────
bin_measures = {'EDSS',  'MSPro'};
bin_thresh   = [3,        1     ];
use_geq      = false;   % y = score > threshold  (matches paper_multi_bi_ridge.m)

cont_measures = {'T25FW', 'x9HPTD', 'x9HPTND', 'SDMTcorrect', 'MSFC_SDMT'};
cont_labels   = {'T25FW', '9HPT-D', '9HPT-ND', 'SDMT', 'MSFC-SDMT'};
cont_transforms = {'log', 'log', 'log', 'raw', 'raw-signed'};

if all_metrics_mode
    qmri_metrics = {'T1', 'MTR', 'FA', 'MD', 'NDI', 'ODI', 'FWF'};
else
    qmri_metrics = {'NDI', 'ODI', 'FWF'};
end
legacy_metrics = {'T1', 'MTR', 'FA', 'MD'};
noddi_metrics = {'NDI', 'ODI', 'FWF'};

% NRF predictor names — the metric is encoded in the GroupTract filename,
% not in the column names. Original tables are .xlsx; NODDI tables are .csv.
% Groups verified from actual file headers.
tract_groups_full = {'Association','Cerebellar','Occipitoparietal','ProjectionBrainstem'};

ico_regions = {'periventricular','juxtacortical','infratentorial','deepwhitematter'};

% NRF predictor list (same for every metric — metric encoded in the file)
nrf_preds_base = build_nrf_preds(tract_groups_full);   % 16 names

demo_cols_bin  = {'Age','GenderNum','DurationOfDisease'};  % matches paper_multi_bi_ridge.m default
demo_cols_cont = {'Age','GenderNum'};                       % matches paper_multi_continuous_ridge_all.m default
n_boot    = p.Results.nBoot;
alpha_ci  = 0.05;
family_size = p.Results.globalFamilySize;

%% ── Load clinical data ────────────────────────────────────────────────────
fprintf('Loading clinical data...\n');
require_file(fullfile(stats_dir, 'clinicalScore.xlsx'));
T_clin = readtable(fullfile(stats_dir, 'clinicalScore.xlsx'));
[T_clin, scoring_audit] = noddi_reproduction_support( ...
    'add_msfc_sdmt', T_clin);

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

%% ── Load classical-region NODDI data (one shared table) ─────────────
% Columns are {region}{metric}; the extraction table keeps the standard AMICO
% free-water name FWF internally, while figures label it ISOVF.
fprintf('Loading classical-region data...\n');
require_file(classical_file);
T_ico_noddi = readtable(classical_file);
validate_classical_table(T_ico_noddi, ico_regions, noddi_metrics);
classical_model_source = classical_file;
if all_metrics_mode
    % A full reproduction run actively regenerates one seven-metric classical
    % table from the source maps. Prefer that table so the models genuinely
    % consume the new MATLAB aggregation. The legacy merge remains as a
    % compatibility fallback for statistics-only reruns made from the original
    % manuscript workbooks plus the historical NODDI CSV.
    legacy_classical_file = fullfile(stats_dir, 'icometrixLesionStats.xlsx');
    require_file(legacy_classical_file);
    T_ico_legacy = readtable(legacy_classical_file);
    validate_legacy_classical_table(T_ico_legacy, ico_regions, legacy_metrics);
    combined_classical_file = fullfile(metrics_dir, ...
        'ClassicalRegionAllMetrics.csv');
    if isfile(combined_classical_file)
        T_ico = readtable(combined_classical_file);
        validate_combined_classical_table(T_ico, ico_regions, ...
            [legacy_metrics noddi_metrics]);
        assert_classical_noddi_identity(T_ico, T_ico_noddi, ...
            ico_regions, noddi_metrics);
        classical_legacy_validation = ...
            validate_combined_classical_legacy_identity( ...
                T_ico, T_ico_legacy, ico_regions, legacy_metrics);
        classical_model_source = combined_classical_file;
    else
        assert(isequal(sort(string(T_ico_legacy.SubjectID)), ...
                       sort(string(T_ico_noddi.SubjectID))), ...
            'Legacy and NODDI classical tables must contain the same subjects.');
        T_ico = innerjoin(T_ico_legacy, T_ico_noddi, 'Keys','SubjectID');
        classical_legacy_validation = struct( ...
            'status','PASS', 'mode','historical table used directly', ...
            'rows_checked',height(T_ico_legacy), ...
            'cells_checked',height(T_ico_legacy) * ...
                numel(ico_regions) * numel(legacy_metrics), ...
            'maximum_absolute_error',0, 'missingness_mismatches',0);
    end
else
    T_ico = T_ico_noddi;
end
T_crf_full = outerjoin(T_clin, T_ico, ...
    'Keys','SubjectID','MergeKeys',true);
T_crf_ms = T_crf_full(~startsWith(T_crf_full.SubjectID, 'sub-C'), :);

%% ── Load NRF data per metric (GroupTract{metric}_All.csv) ────────────────
% Every metric uses the same 16 predictor names.
fprintf('Loading NRF data per metric...\n');
T_nrf_full = struct();   % controls retained  (for continuous)
T_nrf_ms   = struct();   % controls removed   (for binary)
for mi = 1:numel(qmri_metrics)
    metric   = qmri_metrics{mi};
    if ismember(metric, legacy_metrics)
        generated_nrf_file = fullfile(metrics_dir, ...
            sprintf('GroupTract%s_All.xlsx', metric));
        if isfile(generated_nrf_file)
            nrf_file = generated_nrf_file;
        else
            nrf_file = fullfile(stats_dir, ...
                sprintf('GroupTract%s_All.xlsx', metric));
        end
    else
        nrf_file = fullfile(metrics_dir, sprintf('GroupTract%s_All.csv', metric));
    end
    if ~exist(nrf_file, 'file')
        error('NRF file not found: %s', nrf_file);
    end
    T_tract = readtable(nrf_file);
    validate_tract_table(T_tract, nrf_preds_base, metric, ...
        ismember(metric, noddi_metrics));
    T_full  = outerjoin(T_clin, T_tract, 'Keys','SubjectID','MergeKeys',true);
    T_nrf_full.(metric) = T_full;
    T_nrf_ms.(metric)   = T_full(~startsWith(T_full.SubjectID, 'sub-C'), :);
    fprintf('  %s: %d subjects total, %d MS patients\n', metric, ...
        height(T_full), height(T_nrf_ms.(metric)));
end

nrf_preds_base = validate_preds(nrf_preds_base, T_nrf_full.(qmri_metrics{1}), 'NRF-qMRI');

%% ── Output ────────────────────────────────────────────────────────────────
results  = struct();
all_rows = struct([]);   % each metric/outcome → no-demo and +Demo rows
full_figure_rows = struct([]);     % continuous Figure 6/S2 input rows
matched_figure_rows = struct([]);  % exact matched-subject refits

report = struct();
report.schema_version = 4;
report.analysis_scope = analysis_scope;
report.analysis_source = 'MS100Public/rev_qMRI.m';
report.analysis_source_sha256 = ...
    'df4feae5033fbed00b5d749140560be9ddc516333f9552da1700d4abf3bd8c45';
report.source_derived_script = 'helper/matlab/rev_qMRI_noddi.m';
if all_metrics_mode
    report.entry_point = 'code/run_all_qmri_manuscript_models.m';
    report.adaptation = [ ...
        'Unified active refit of T1/MTR/FA/MD/NDI/ODI/FWF with raw ' ...
        'SDMT/signed raw MSFC-SDMT outcomes; original ridge mechanics retained.'];
else
    report.entry_point = 'helper/matlab/rev_qMRI_noddi.m';
    report.adaptation = [ ...
        'NODDI-only compatibility analysis with raw SDMT/signed raw ' ...
        'MSFC-SDMT outcomes; original ridge mechanics retained.'];
end
report.edss_rule = 'EDSS > 3.0 (strict; current manuscript rule)';
report.mspro_rule = 'MSPro > 1 (equivalent to >=2 for integer MSPro)';
report.noddi_case_rule = [ ...
    'binary joint tract+classical complete cases; continuous independent ' ...
    'available cases; separate joint-case refits for figure markers'];
report.classical_input = noddi_reproduction_support( ...
    'classical_input_provenance', classical_file, ...
    fullfile(stats_dir,'icometrixLesionStats.xlsx'));
classical_ids = string(T_ico_noddi.SubjectID);
clinical_ids = string(T_clin.SubjectID);
report.classical_input.clinical_overlap = numel(intersect(classical_ids,clinical_ids));
report.classical_input.classical_only_subjects = numel(setdiff(classical_ids,clinical_ids));
report.classical_input.clinical_only_subjects = numel(setdiff(clinical_ids,classical_ids));
report.classical_input.schema_validated = true;
report.classical_input.display_alias = 'FWF is displayed as ISOVF';
report.classical_model_source = classical_model_source;
if all_metrics_mode
    report.classical_legacy_validation = classical_legacy_validation;
end
report.msfc = scoring_audit;
report.inputs = noddi_reproduction_support( ...
    'input_provenance', stats_dir, metrics_dir);
input_matches = [report.inputs.matches_completed_run];
if ~all(input_matches)
    mismatched_names = {report.inputs(~input_matches).name};
    error('Required aggregate input fingerprint mismatch: %s', ...
        strjoin(mismatched_names, ', '));
end
report.all_input_fingerprints_match_completed_run = true;
report.cohort_validation = noddi_reproduction_support( ...
    'validate_input_cohorts', T_clin, metrics_dir);

%% ── Legacy FA/MD equivalence gate ────────────────────────────────────────
fprintf('\nVALIDATION: refitting tract FA/MD with original case rules...\n');
legacy_rows = run_legacy_validation(T_clin, stats_dir, metrics_dir, nrf_preds_base, ...
    ico_regions, bin_measures, bin_thresh, cont_measures(1:3), ...
    cont_labels(1:3), demo_cols_bin, demo_cols_cont, n_boot, alpha_ci);
validation = noddi_reproduction_support('compare_legacy_rows', legacy_rows);
if validation.all_match
    validation.status = 'PASS';
else
    validation.recovered_mask_migration = noddi_reproduction_support( ...
        'validate_recovered_mask_migration', validation);
    validation.status = 'PASS_EXPECTED_RECOVERED_MASK_MIGRATION';
end
validation_path = fullfile(output_dir, 'fa_md_equivalence_validation.json');
noddi_reproduction_support('write_json', validation_path, validation);
report.legacy_validation = validation;
fprintf('PASS: historical FA/MD gate status %s (%d/%d exact rows).\n', ...
    validation.status, validation.rows_matched, validation.rows_checked);

%% ════════════════════════════════════════════════════════════════════════
%  BINARY OUTCOMES
%  Joint complete cases: subjects need finite data in BOTH tract and
%  classical predictor columns. Exact zero is a valid NODDI estimate and is
%  retained. Both ridge fits precede either AUC evaluation so the second fit
%  inherits the same RNG state as the original tract/classical pair.
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
        exclude_predictor_zeros = ismember(metric, legacy_metrics);
        nrf_preds = validate_preds(nrf_preds_base, T_nrf_ms.(metric), sprintf('NRF-%s', metric));
        crf_preds = validate_preds(build_crf_preds(ico_regions, metric), ...
                                   T_crf_ms, sprintf('CRF-%s', metric));
        assert(numel(nrf_preds) == 16 && numel(crf_preds) == 4, ...
            'Expected 16 tract and four classical predictors for %s.', metric);

        name_nrf = sprintf('Tract-based %s', metric);
        name_crf = sprintf('Classical-region %s', metric);

        % ── No-demo: no demographic filtering (demographics={} in original)
        [Xn_nd, Xc_nd, y_nd, ~, n_nd] = get_joint_cases( ...
            T_nrf_ms.(metric), T_crf_ms, nrf_preds, crf_preds, ...
            measure, {}, true, threshold, use_geq, ...
            exclude_predictor_zeros, 'binary');

        lambda_grid_bin = logspace(-6, 6, 50);   % 50-point grid, matches paper_multi_bi_ridge.m

        % ── No-demo pair ──────────────────────────────────────────────────
        if n_nd >= 20
            rng(42);
            [LP_nrf_nd, lambda_nd] = fit_ridge_lp_bin(zscore(Xn_nd), y_nd, lambda_grid_bin, ...
                                         sprintf('%s|%s|nd', measure, name_nrf));
            [LP_crf_nd, lambda_crf_nd] = fit_ridge_lp_bin(zscore(Xc_nd), y_nd, lambda_grid_bin, ...
                                         sprintf('%s|%s|nd', measure, name_crf));

            res = eval_binary(LP_nrf_nd, y_nd, [], n_nd, n_boot, alpha_ci, ...
                              sprintf('%s|%s|nodemo', measure, name_nrf));
            results.(sprintf('%s_NRF_%s_nd', measure, metric)) = res;
            all_rows = append_canonical_row(all_rows, measure, measure, ...
                metric, 'Tract-based', false, {}, ...
                'joint tract+classical complete cases', res, lambda_nd, ...
                true, y_nd, family_size);

            res = eval_binary(LP_crf_nd, y_nd, [], n_nd, n_boot, alpha_ci, ...
                              sprintf('%s|%s|nodemo', measure, name_crf));
            results.(sprintf('%s_CRF_%s_nd', measure, metric)) = res;
            all_rows = append_canonical_row(all_rows, measure, measure, ...
                metric, 'Classical-region', false, {}, ...
                'joint tract+classical complete cases', res, lambda_crf_nd, ...
                true, y_nd, family_size);
        else
            warn_skip(measure, sprintf('%s/%s no-demo', name_nrf, name_crf), n_nd);
        end

        % ── With-demo pair ────────────────────────────────────────────────
        [Xn_wd, Xc_wd, y_wd, Demo_wd, n_wd] = get_joint_cases( ...
            T_nrf_ms.(metric), T_crf_ms, nrf_preds, crf_preds, ...
            measure, demo_cols_bin, true, threshold, use_geq, ...
            exclude_predictor_zeros, 'binary');

        if n_wd >= 20
            rng(42);
            [LP_nrf_wd, lambda_wd] = fit_ridge_lp_bin(zscore(Xn_wd), y_wd, lambda_grid_bin, ...
                                         sprintf('%s|%s|wd', measure, name_nrf));
            [LP_crf_wd, lambda_crf_wd] = fit_ridge_lp_bin(zscore(Xc_wd), y_wd, lambda_grid_bin, ...
                                         sprintf('%s|%s|wd', measure, name_crf));

            res = eval_binary(LP_nrf_wd, y_wd, Demo_wd, n_wd, n_boot, alpha_ci, ...
                              sprintf('%s|%s|demo', measure, name_nrf));
            results.(sprintf('%s_NRF_%s_wd', measure, metric)) = res;
            all_rows = append_canonical_row(all_rows, measure, measure, ...
                metric, 'Tract-based', true, demo_cols_bin, ...
                'joint tract+classical complete cases', res, lambda_wd, ...
                true, y_wd, family_size);

            res = eval_binary(LP_crf_wd, y_wd, Demo_wd, n_wd, n_boot, alpha_ci, ...
                              sprintf('%s|%s|demo', measure, name_crf));
            results.(sprintf('%s_CRF_%s_wd', measure, metric)) = res;
            all_rows = append_canonical_row(all_rows, measure, measure, ...
                metric, 'Classical-region', true, demo_cols_bin, ...
                'joint tract+classical complete cases', res, lambda_crf_wd, ...
                true, y_wd, family_size);
        else
            warn_skip(measure, sprintf('%s/%s +demo', name_nrf, name_crf), n_wd);
        end
    end
end

%% ════════════════════════════════════════════════════════════════════════
%  CONTINUOUS OUTCOMES
%  Independent tract and classical available cases, matching the original
%  model1/model2 behavior. Separately refit both branches on their exact joint
%  MS-only cases for Figure 6/S2 matched-subject markers.
%% ════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n  CONTINUOUS OUTCOMES\n%s\n', repmat('=',1,65), repmat('=',1,65));

for ci = 1:numel(cont_measures)
    measure = cont_measures{ci};
    clabel  = cont_labels{ci};
    transform = cont_transforms{ci};
    fprintf('\n>>> %s\n', clabel);

    for mi = 1:numel(qmri_metrics)
        metric    = qmri_metrics{mi};
        exclude_predictor_zeros = ismember(metric, legacy_metrics);
        nrf_preds = validate_preds(nrf_preds_base, T_nrf_full.(metric), sprintf('NRF-%s', metric));
        crf_preds = validate_preds(build_crf_preds(ico_regions, metric), ...
                                   T_crf_full, sprintf('CRF-%s', metric));
        assert(numel(nrf_preds) == 16 && numel(crf_preds) == 4, ...
            'Expected 16 tract and four classical predictors for %s.', metric);
        name_nrf = sprintf('Tract-based %s', metric);
        name_crf = sprintf('Classical-region %s', metric);

        % Load all independent case sets upfront (data loading is deterministic).
        [Xn_nd, yn_nd, ~,    nn_nd] = get_complete_cases( ...
            T_nrf_full.(metric), nrf_preds, measure, {}, false, [], [], ...
            exclude_predictor_zeros, transform);
        [Xn_wd, yn_wd, Dn_wd, nn_wd] = get_complete_cases( ...
            T_nrf_full.(metric), nrf_preds, measure, demo_cols_cont, ...
            false, [], [], exclude_predictor_zeros, transform);
        [Xc_nd, yc_nd, ~,    nc_nd] = get_complete_cases( ...
            T_crf_full, crf_preds, measure, {}, false, [], [], ...
            exclude_predictor_zeros, transform);
        [Xc_wd, yc_wd, Dc_wd, nc_wd] = get_complete_cases( ...
            T_crf_full, crf_preds, measure, demo_cols_cont, ...
            false, [], [], exclude_predictor_zeros, transform);

        rng(42);
        if nn_nd >= 20
            res_n_nd = run_model_cont(Xn_nd, yn_nd, [], nn_nd, ...
                sprintf('%s|%s|nd', clabel, name_nrf));
            results.(sprintf('%s_NRF_%s_nd', measure, metric)) = res_n_nd;
            all_rows = append_canonical_row(all_rows, clabel, measure, ...
                metric, 'Tract-based', false, {}, ...
                'independent tract available cases', res_n_nd, ...
                res_n_nd.lambda, false, [], family_size);
        else
            warn_skip(clabel, sprintf('%s no-demo', name_nrf), nn_nd);
        end
        if nc_nd >= 20
            res_c_nd = run_model_cont(Xc_nd, yc_nd, [], nc_nd, ...
                sprintf('%s|%s|nd', clabel, name_crf));
            results.(sprintf('%s_CRF_%s_nd', measure, metric)) = res_c_nd;
            all_rows = append_canonical_row(all_rows, clabel, measure, ...
                metric, 'Classical-region', false, {}, ...
                'independent classical available cases', res_c_nd, ...
                res_c_nd.lambda, false, [], family_size);
        else
            warn_skip(clabel, sprintf('%s no-demo', name_crf), nc_nd);
        end

        rng(42);
        if nn_wd >= 20
            res_n_wd = run_model_cont(Xn_wd, yn_wd, Dn_wd, nn_wd, ...
                sprintf('%s|%s|wd', clabel, name_nrf));
            results.(sprintf('%s_NRF_%s_wd', measure, metric)) = res_n_wd;
            all_rows = append_canonical_row(all_rows, clabel, measure, ...
                metric, 'Tract-based', true, demo_cols_cont, ...
                'independent tract available cases', res_n_wd, ...
                res_n_wd.lambda, false, [], family_size);
        else
            warn_skip(clabel, sprintf('%s +demo', name_nrf), nn_wd);
        end
        if nc_wd >= 20
            res_c_wd = run_model_cont(Xc_wd, yc_wd, Dc_wd, nc_wd, ...
                sprintf('%s|%s|wd', clabel, name_crf));
            results.(sprintf('%s_CRF_%s_wd', measure, metric)) = res_c_wd;
            all_rows = append_canonical_row(all_rows, clabel, measure, ...
                metric, 'Classical-region', true, demo_cols_cont, ...
                'independent classical available cases', res_c_wd, ...
                res_c_wd.lambda, false, [], family_size);
        else
            warn_skip(clabel, sprintf('%s +demo', name_crf), nc_wd);
        end

        assert(exist('res_n_nd','var') && exist('res_c_nd','var') && ...
            exist('res_n_wd','var') && exist('res_c_wd','var'), ...
            'Full continuous model pair was not produced for %s/%s.', clabel, metric);
        full_figure_rows = append_figure_pair(full_figure_rows, measure, metric, ...
            res_n_nd, res_n_wd, res_c_nd, res_c_wd);

        % Exact matched-subject refits use MS-only rows and joint completeness.
        [Xnm_nd, Xcm_nd, ym_nd, ~, nm_nd] = get_joint_cases( ...
            T_nrf_ms.(metric), T_crf_ms, nrf_preds, crf_preds, ...
            measure, {}, false, [], [], exclude_predictor_zeros, transform);
        [Xnm_wd, Xcm_wd, ym_wd, Dm_wd, nm_wd] = get_joint_cases( ...
            T_nrf_ms.(metric), T_crf_ms, nrf_preds, crf_preds, ...
            measure, demo_cols_cont, false, [], [], ...
            exclude_predictor_zeros, transform);
        assert(nm_nd >= 20 && nm_wd >= 20, ...
            'Too few exact matched cases for %s/%s.', clabel, metric);

        rng(42);
        res_nm_nd = run_model_cont(Xnm_nd, ym_nd, [], nm_nd, ...
            sprintf('%s|%s|matched-nd', clabel, name_nrf));
        res_cm_nd = run_model_cont(Xcm_nd, ym_nd, [], nm_nd, ...
            sprintf('%s|%s|matched-nd', clabel, name_crf));
        rng(42);
        res_nm_wd = run_model_cont(Xnm_wd, ym_wd, Dm_wd, nm_wd, ...
            sprintf('%s|%s|matched-wd', clabel, name_nrf));
        res_cm_wd = run_model_cont(Xcm_wd, ym_wd, Dm_wd, nm_wd, ...
            sprintf('%s|%s|matched-wd', clabel, name_crf));
        matched_figure_rows = append_figure_pair(matched_figure_rows, ...
            measure, metric, res_nm_nd, res_nm_wd, res_cm_nd, res_cm_wd);

        clear res_n_nd res_c_nd res_n_wd res_c_wd
    end
end

%% ── Print & save ──────────────────────────────────────────────────────────
model_table = struct2table(all_rows);
expected_rows = numel(qmri_metrics) * 7 * 2 * 2;
expected_adjusted = expected_rows / 2;
expected_per_framework = expected_rows / 2;
assert(height(model_table) == expected_rows, ...
    'Expected %d model rows; found %d.', expected_rows, height(model_table));
assert(sum(model_table.Adjusted) == expected_adjusted, ...
    'Expected %d adjusted rows.', expected_adjusted);
assert(nnz(strcmp(model_table.Model,'Tract-based')) == expected_per_framework && ...
       nnz(strcmp(model_table.Model,'Classical-region')) == expected_per_framework, ...
    'Expected %d rows for each model framework.', expected_per_framework);

if all_metrics_mode
    model_filename = 'all_qmri_manuscript_models.csv';
    if strcmp(reference_validation_mode, 'frozen')
        report.unified_reference_validation = ...
            validate_unified_reference_results(model_table);
    else
        report.unified_reference_validation = struct( ...
            'status','DEFERRED_ONE_TIME_REFERENCE_MIGRATION', ...
            'reason',[ ...
                'Generate the recovered-132 candidate once, freeze its reviewed ' ...
                'references, then rerun in frozen mode.']);
    end
else
    model_filename = 'noddi_manuscript_aligned_models.csv';
end
out_file = fullfile(output_dir, model_filename);
writetable(model_table, out_file);
fprintf('\nResults saved to: %s\n', out_file);

if all_metrics_mode
    noddi_compatibility = model_table( ...
        ismember(string(model_table.Metric), ["NDI","ODI","FWF"]), :);
    assert(height(noddi_compatibility) == 84, ...
        'Expected 84 rows in the NODDI compatibility subset.');
    noddi_compatibility_filename = 'noddi_manuscript_aligned_models.csv';
    writetable(noddi_compatibility, ...
        fullfile(output_dir, noddi_compatibility_filename));
end

[figure_full_path, figure_matched_path] = write_scope_figure_model_tables( ...
    full_figure_rows, matched_figure_rows, output_dir, all_metrics_mode);

report.results_file = model_filename;
report.model_rows = height(model_table);
report.metrics = qmri_metrics;
report.primary_demo_rows = sum(model_table.Adjusted);
if ~all_metrics_mode
    % Retain fields consumed by the frozen NODDI-only reproduction report.
    report.noddi_results_file = model_filename;
    report.noddi_rows = height(model_table);
end
report.global_family_size = family_size;
report.global_family_definition = family_definition(family_size);
report.no_demo_multiplicity = [ ...
    'No-demo rows retain raw p only. The expanded 98-test manuscript family ' ...
    'uses the +Demo rows; it does not include a second no-demo family.'];
if all_metrics_mode
    report_filename = 'all_qmri_manuscript_analysis_report.json';
else
    report_filename = 'noddi_manuscript_analysis_report.json';
end
report.output_files = struct( ...
    'models', model_filename, ...
    'figure6_full_models', get_filename(figure_full_path), ...
    'figure6_matched_models', get_filename(figure_matched_path), ...
    'validation', 'fa_md_equivalence_validation.json', ...
    'report', report_filename);
if all_metrics_mode
    report.output_files.noddi_compatibility_models = ...
        noddi_compatibility_filename;
end
noddi_reproduction_support('write_json', ...
    fullfile(output_dir, report.output_files.report), report);
end  % ── main ──────────────────────────────────────────────────────────────


%% ══════════════════════════════════════════════════════════════════════════
%  ANALYSIS FUNCTIONS
%% ══════════════════════════════════════════════════════════════════════════

function [LP, lambda] = fit_ridge_lp_bin(X_z, y, lambda_grid, tag)
% Phase 1: fit ridge logistic and return the linear predictor LP.
%
% rng is NOT reset here. The caller sets rng(42) immediately before this
% tract-model fit, retaining the deterministic behavior of rev_qMRI.m.

    cv    = fitclinear(X_z, y, 'Learner','logistic', ...
                'Regularization','ridge', 'Lambda',lambda_grid, 'KFold',10);
    [~, bi] = min(kfoldLoss(cv));
    lambda = lambda_grid(bi);
    mdl_r   = fitclinear(X_z, y, 'Learner','logistic', ...
                'Regularization','ridge', 'Lambda',lambda);
    LP = X_z * mdl_r.Beta + mdl_r.Bias;
    fprintf('  [ridge] %s  lambda=%.2e\n', tag, lambda);
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
% y is already transformed as specified by the outcome: log for T25FW/9HPT,
% raw for SDMT, and signed raw for MSFC-SDMT.
%
% rng is NOT reset here; the caller resets rng(42) immediately before the
% tract fit, retaining the deterministic behavior of rev_qMRI.m.

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

    res = struct('n',n, 'p',p, 'AIC',AIC, 'R2',R2, ...
        'lambda',lambda_grid(bi));
end


%% ══════════════════════════════════════════════════════════════════════════
%  ROW BUILDERS  (one row per call)
%% ══════════════════════════════════════════════════════════════════════════

function [X, y, Demo, n] = get_complete_cases(T, pred_names, measure, ...
                                               demo_cols, is_binary, threshold, use_geq, ...
                                               excl_pred_zeros, transform)
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
        valid = pred_ok & isfinite(y_raw) & demo_ok;
        if strcmp(transform, 'log')
            valid = valid & y_raw > 0;
        end
    end

    X    = X_raw(valid, :);
    Demo = Dmat(valid, :);   % Nx0 when demo_cols={}; isempty(Demo) → true
    n    = sum(valid);
    y_v  = y_raw(valid);
    if is_binary
        if use_geq, y = double(y_v >= threshold); else, y = double(y_v > threshold); end
    elseif strcmp(transform, 'log')
        y = log(y_v);
    else
        y = double(y_v);
    end
end


function [X_nrf, X_crf, y, Demo, n] = get_joint_cases(T_nrf, T_crf, ...
        nrf_preds, crf_preds, measure, demo_cols, is_binary, threshold, use_geq, ...
        excl_pred_zeros, transform)
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
        valid = pred_ok & isfinite(y_raw) & demo_ok;
        if strcmp(transform, 'log')
            valid = valid & y_raw > 0;
        end
    end

    X_nrf = X_nrf_raw(valid, :);
    X_crf = X_crf_raw(valid, :);
    Demo  = Dmat(valid, :);   % Nx0 when demo_cols={}
    n     = sum(valid);
    y_v   = y_raw(valid);
    if is_binary
        if use_geq, y = double(y_v >= threshold); else, y = double(y_v > threshold); end
    elseif strcmp(transform, 'log')
        y = log(y_v);
    else
        y = double(y_v);
    end
end


%% ══════════════════════════════════════════════════════════════════════════
%  PREDICTOR BUILDERS
%% ══════════════════════════════════════════════════════════════════════════

function preds = build_nrf_preds(tract_groups_full)
% 16 predictors: 4 groups × L/R × {Tail, NAWM}.
% Column names verified from the GroupTract{metric}_All tables.
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

function validate_classical_table(T, ico_regions, qmri_metrics)
% Require the exact extraction schema so a renamed/misordered private table
% cannot silently change the four-predictor classical models.
    expected = {'SubjectID'};
    for ri = 1:numel(ico_regions)
        for mi = 1:numel(qmri_metrics)
            expected{end+1} = [ico_regions{ri} qmri_metrics{mi}]; %#ok<AGROW>
        end
    end
    assert(isequal(T.Properties.VariableNames, expected), ...
        'Classical NODDI table schema/order does not match the required 13 columns.');
    subjects = string(T.SubjectID);
    assert(all(~ismissing(subjects) & strlength(strtrim(subjects)) > 0), ...
        'Classical NODDI table contains a missing SubjectID.');
    assert(numel(unique(subjects)) == height(T), ...
        'Classical NODDI table must contain one row per SubjectID.');
    assert(height(T) == 87, ...
        'Expected the 87-subject legacy classical-region cohort; found %d.', height(T));
    assert(~any(startsWith(subjects,'sub-C')), ...
        'Classical NODDI table must contain MS participants only.');
    values = table2array(T(:,2:end));
    assert(isnumeric(values) && all(isfinite(values) | isnan(values),'all'), ...
        'Classical NODDI values must be numeric finite values or NaN.');
    observed = values(isfinite(values));
    assert(all(observed >= 0 & observed <= 1), ...
        'Classical NODDI values must lie in [0,1].');
end

function validate_legacy_classical_table(T, ico_regions, legacy_metrics)
% The historical table also contains RD/AD. Require every column used by
% this unified refit without changing or reordering the legacy workbook.
    expected = {'SubjectID'};
    for ri = 1:numel(ico_regions)
        for mi = 1:numel(legacy_metrics)
            expected{end+1} = [ico_regions{ri} legacy_metrics{mi}]; %#ok<AGROW>
        end
    end
    missing = expected(~ismember(expected, T.Properties.VariableNames));
    assert(isempty(missing), ...
        'Legacy classical table is missing columns: %s', strjoin(missing, ', '));
    subjects = string(T.SubjectID);
    assert(height(T) == 87 && numel(unique(subjects)) == 87, ...
        'Legacy classical table must contain 87 unique subjects.');
    assert(all(~ismissing(subjects) & strlength(strtrim(subjects)) > 0), ...
        'Legacy classical table contains a missing SubjectID.');
    assert(~any(startsWith(subjects,'sub-C')), ...
        'Legacy classical table must contain MS participants only.');
    values = table2array(T(:,expected(2:end)));
    assert(isnumeric(values) && all(isfinite(values) | isnan(values),'all'), ...
        'Legacy classical predictors must be numeric finite values or NaN.');
end


function validate_combined_classical_table(T, ico_regions, metrics)
% Require the exact region-major schema emitted by generate_all_metric_tables.
expected = {'SubjectID'};
for ri = 1:numel(ico_regions)
    for mi = 1:numel(metrics)
        expected{end+1} = [ico_regions{ri} metrics{mi}]; %#ok<AGROW>
    end
end
assert(isequal(T.Properties.VariableNames, expected), ...
    'Combined classical table schema/order is not SubjectID plus 4 x 7 metrics.');
subjects = string(T.SubjectID);
assert(height(T) == 87 && numel(unique(subjects)) == 87, ...
    'Combined classical table must contain 87 unique subjects.');
assert(all(~ismissing(subjects) & strlength(strtrim(subjects)) > 0), ...
    'Combined classical table contains a missing SubjectID.');
assert(~any(startsWith(subjects,'sub-C')), ...
    'Combined classical table must contain MS participants only.');
values = table2array(T(:,2:end));
assert(isnumeric(values) && all(isfinite(values) | isnan(values),'all'), ...
    'Combined classical predictors must be numeric finite values or NaN.');
noddiNames = expected(contains(expected, {'NDI','ODI','FWF'}));
noddiValues = table2array(T(:,noddiNames));
observed = noddiValues(isfinite(noddiValues));
assert(all(observed >= 0 & observed <= 1), ...
    'Combined classical NODDI values must lie in [0,1].');
end


function assert_classical_noddi_identity(combined, noddi, ico_regions, noddi_metrics)
% The combined table is written in the same MATLAB extraction call as the
% compatibility NODDI table. Require exact agreement (including missingness)
% before allowing it to become the modelling input.
combined_ids = string(combined.SubjectID);
noddi_ids = string(noddi.SubjectID);
[present, order] = ismember(combined_ids, noddi_ids);
assert(all(present) && height(combined) == height(noddi), ...
    'Combined and NODDI classical tables must contain identical subjects.');
for ri = 1:numel(ico_regions)
    for mi = 1:numel(noddi_metrics)
        name = [ico_regions{ri} noddi_metrics{mi}];
        noddi_values = noddi.(name);
        assert(isequaln(combined.(name), noddi_values(order)), ...
            'Combined classical table differs from the NODDI table in %s.', name);
    end
end
end


function audit = validate_combined_classical_legacy_identity( ...
        combined, legacy, ico_regions, legacy_metrics)
% Bind every legacy-valued cell in the active seven-metric table to the
% original 87-subject table before the combined table enters a model.
assert(height(combined) == height(legacy) && ...
    isequal(string(combined.SubjectID), string(legacy.SubjectID)), ...
    'Combined and historical classical SubjectID rows/order differ.');
maximum_error = 0;
cells_checked = 0;
missingness_mismatches = 0;
tolerance = 1e-12;
for ri = 1:numel(ico_regions)
    for mi = 1:numel(legacy_metrics)
        name = [ico_regions{ri} legacy_metrics{mi}];
        observed = double(combined.(name));
        expected = double(legacy.(name));
        mismatch = xor(isnan(observed), isnan(expected));
        missingness_mismatches = missingness_mismatches + nnz(mismatch);
        finite = isfinite(observed) & isfinite(expected);
        difference = abs(observed(finite) - expected(finite));
        if ~isempty(difference)
            maximum_error = max(maximum_error, max(difference));
        end
        cells_checked = cells_checked + numel(observed);
        assert(~any(mismatch) && all(difference <= tolerance), ...
            'Combined classical legacy values differ in %s.', name);
    end
end
audit = struct('status','PASS', ...
    'mode','direct legacy-column comparison', ...
    'rows_checked',height(combined), 'cells_checked',cells_checked, ...
    'maximum_absolute_error',maximum_error, ...
    'missingness_mismatches',missingness_mismatches, ...
    'absolute_tolerance',tolerance);
end


function validate_tract_table(T, predictor_names, metric, is_noddi)
% Validate the analyzed 16 columns and the cohort before any clinical join.
    required = [{'SubjectID'}, predictor_names];
    missing = required(~ismember(required, T.Properties.VariableNames));
    assert(isempty(missing), '%s tract table is missing columns: %s', ...
        metric, strjoin(missing, ', '));
    subjects = string(T.SubjectID);
    assert(all(~ismissing(subjects) & strlength(strtrim(subjects)) > 0), ...
        '%s tract table contains a missing SubjectID.', metric);
    assert(numel(unique(subjects)) == height(T), ...
        '%s tract table must contain one row per SubjectID.', metric);
    expected_rows = 132;
    assert(height(T) == expected_rows, ...
        '%s tract table must contain %d subjects; found %d.', ...
        metric, expected_rows, height(T));
    values = table2array(T(:,predictor_names));
    assert(isnumeric(values) && all(isfinite(values) | isnan(values),'all'), ...
        '%s tract predictors must be numeric finite values or NaN.', metric);
    if is_noddi
        observed = values(isfinite(values));
        assert(all(observed >= 0 & observed <= 1), ...
            '%s tract predictors must lie in [0,1].', metric);
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

function warn_skip(outcome, model, n)
    warning('rev_qMRI: only %d cases for %s — %s. Skipping.', n, outcome, model);
end

function s = pval_str(p)
    if p < 0.001, s = '<0.001'; else, s = sprintf('%.4f', p); end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end


%% ══════════════════════════════════════════════════════════════════════════
%  REVIEWER-PACKAGE ROWS AND LEGACY EQUIVALENCE
%% ══════════════════════════════════════════════════════════════════════════

function rows = append_canonical_row(rows, outcome, source_field, metric, ...
        model_name, adjusted, demo_cols, case_policy, res, lambda, ...
        is_binary, y, family_size)
% Add the canonical CSV fields without changing the fitted-model mechanics.
    row = struct();
    row.Outcome = outcome;
    row.SourceField = source_field;
    row.Metric = metric;
    row.Model = model_name;
    row.Adjusted = adjusted;
    row.Demographics = strjoin(demo_cols, '+');
    if isempty(row.Demographics), row.Demographics = 'none'; end
    row.CasePolicy = case_policy;
    row.N = res.n;
    if is_binary
        row.PositiveN = sum(y);
        row.PerformanceName = 'apparent AUC';
        row.Performance = res.AUC_mean;
        row.CI_lower = res.AUC_CI(1);
        row.CI_upper = res.AUC_CI(2);
    else
        row.PositiveN = NaN;
        row.PerformanceName = 'apparent R2';
        row.Performance = res.R2;
        row.CI_lower = NaN;
        row.CI_upper = NaN;
    end
    row.Lambda = lambda;
    row.p_raw = res.p;
    if adjusted && isfinite(family_size)
        row.p_adj = min(res.p * family_size, 1);
        row.MultiplicityFamilyN = family_size;
        row.InferenceRole = 'primary +Demo qMRI family';
    else
        row.p_adj = NaN;
        row.MultiplicityFamilyN = 0;
        row.InferenceRole = 'no-demo performance; outside primary p family';
    end
    row.AIC = res.AIC;
    rows = append_struct(rows, row);
end


function rows = append_figure_pair(rows, source_field, metric, ...
        tract_no_demo, tract_demo, classic_no_demo, classic_demo)
% Figure 6/S2 stores one N plus no-demo and +Demo R2 in each model row.
% Refuse to collapse different demographic case sets into that one-N schema.
    assert(tract_no_demo.n == tract_demo.n, ...
        'Tract no-demo/+Demo N differs for %s/%s figure row.', source_field, metric);
    assert(classic_no_demo.n == classic_demo.n, ...
        'Classical no-demo/+Demo N differs for %s/%s figure row.', source_field, metric);
    display_metric = metric;
    if strcmp(metric,'FWF'), display_metric = 'ISOVF'; end
    key = sprintf('%s-%s', source_field, display_metric);

    row = struct('ClinicalMetric',key, 'Model','tract-based', ...
        'N',tract_no_demo.n, 'R2',tract_no_demo.R2, ...
        'R2_Demo',tract_demo.R2);
    rows = append_struct(rows,row);
    row = struct('ClinicalMetric',key, 'Model','classic', ...
        'N',classic_no_demo.n, 'R2',classic_no_demo.R2, ...
        'R2_Demo',classic_demo.R2);
    rows = append_struct(rows,row);
end


function validation = validate_unified_reference_results(actual)
% Fail closed unless active unified refits reproduce the scientific content
% of the retained consolidated baseline. Upstream table gates bind full values;
% this gate therefore compares exact identities/cohorts/multiplicity and the
% precision actually reported in the manuscript. Lambda is diagnostic only
% because tied cross-validation minima can select adjacent penalties without
% changing the fitted result at reported precision.
    matlab_dir = fileparts(mfilename('fullpath'));
    helper_dir = fileparts(matlab_dir);
    package_root = fileparts(helper_dir);
    absolute_tolerance = 1e-12;

    baseline_path = fullfile(package_root, 'results', 'model_results', ...
        'all_qmri_manuscript_models.csv');
    require_file(baseline_path);
    baseline = readtable(baseline_path, 'TextType','string');
    assert(height(baseline) == 196, ...
        'Expected 196 rows in the consolidated qMRI reference baseline.');
    assert(isequal(actual.Properties.VariableNames, ...
                   baseline.Properties.VariableNames), ...
        'Active and consolidated reference schemas differ.');

    baseline_adjusted = logical_reference_column(baseline.Adjusted);
    legacy = baseline(baseline_adjusted & ...
        ismember(string(baseline.Metric), ["T1","MTR","FA","MD"]), :);
    legacy.p_adj_m56 = min(legacy.p_raw .* 56, 1);
    legacy.p_adj_m98 = legacy.p_adj;
    legacy.significant_m56 = legacy.p_adj_m56 < 0.05;
    legacy.significant_m98 = legacy.p_adj_m98 < 0.05;
    assert(height(legacy) == 56, ...
        ['Expected 56 pre-NODDI adjusted qMRI reference rows ' ...
         '(four legacy metrics across all seven outcomes).']);
    legacy_actual = actual(actual.Adjusted & ...
        ismember(string(actual.Metric), ["T1","MTR","FA","MD"]), :);
    assert(height(legacy_actual) == 56, ...
        'Expected 56 actively refitted legacy-metric adjusted qMRI rows.');
    legacy_order = match_reference_rows(legacy_actual, legacy, false);
    p_difference = abs(legacy_actual.p_raw(legacy_order) - legacy.p_raw);
    p98_difference = abs(legacy_actual.p_adj(legacy_order) - legacy.p_adj_m98);
    legacy_p_display = p_display_matches( ...
        legacy_actual.p_raw(legacy_order), legacy.p_raw);
    legacy_p98_display = p_display_matches( ...
        legacy_actual.p_adj(legacy_order), legacy.p_adj_m98);
    legacy_significance = legacy_actual.p_adj(legacy_order) < 0.05;
    legacy_m56_significance = min(legacy_actual.p_raw(legacy_order) .* 56, 1) < 0.05;
    active_p98 = min(legacy_actual.p_raw(legacy_order) .* 98, 1);
    assert(all(abs(active_p98 - legacy_actual.p_adj(legacy_order)) <= ...
        absolute_tolerance), ...
        'Legacy-metric adjusted rows do not contain internally valid m=98 p values.');
    assert(all(legacy_p_display) && all(legacy_p98_display), ...
        ['Pre-NODDI adjusted p values changed at manuscript precision ' ...
         '(raw mismatches %d; adjusted mismatches %d).'], ...
        nnz(~legacy_p_display), nnz(~legacy_p98_display));
    reference_m98_significance = ...
        logical_reference_column(legacy.significant_m98);
    reference_m56_significance = ...
        logical_reference_column(legacy.significant_m56);
    assert(isequal(legacy_significance, reference_m98_significance) && ...
           isequal(legacy_m56_significance, reference_m56_significance), ...
        'Pre-NODDI adjusted significance decisions changed.');

    noddi = baseline(ismember(string(baseline.Metric), ...
        ["NDI","ODI","FWF"]), :);
    noddi_actual = actual(ismember(string(actual.Metric), ...
        ["NDI","ODI","FWF"]), :);
    assert(height(noddi) == 84 && height(noddi_actual) == 84, ...
        'Expected 84 NODDI rows in both active and frozen tables.');
    assert(isequal(noddi_actual.Properties.VariableNames, ...
                   noddi.Properties.VariableNames), ...
        'Active and frozen NODDI schemas differ.');
    noddi_order = match_reference_rows(noddi_actual, noddi, true);

    text_columns = {'Demographics','CasePolicy','PerformanceName','InferenceRole'};
    for ci = 1:numel(text_columns)
        name = text_columns{ci};
        assert(all(string(noddi_actual.(name)(noddi_order)) == ...
                   string(noddi.(name))), ...
            'NODDI reference mismatch in %s.', name);
    end
    exact_numeric_columns = {'N','PositiveN','MultiplicityFamilyN'};
    max_noddi_difference = 0;
    for ci = 1:numel(exact_numeric_columns)
        name = exact_numeric_columns{ci};
        observed = double(noddi_actual.(name)(noddi_order));
        expected = double(noddi.(name));
        assert(isequaln(observed, expected), ...
            'NODDI exact reference mismatch in %s.', name);
    end

    display_contract = { ...
        'Performance', 3; 'CI_lower', 3; 'CI_upper', 3; 'AIC', 2};
    display_mismatches = 0;
    for ci = 1:size(display_contract,1)
        name = display_contract{ci,1};
        decimals = display_contract{ci,2};
        observed = double(noddi_actual.(name)(noddi_order));
        expected = double(noddi.(name));
        assert(isequal(isnan(observed), isnan(expected)), ...
            'NODDI reference NaN pattern mismatch in %s.', name);
        finite = isfinite(observed) & isfinite(expected);
        difference = abs(observed(finite) - expected(finite));
        if ~isempty(difference)
            max_noddi_difference = max(max_noddi_difference, max(difference));
        end
        matches = numeric_display_matches(observed, expected, decimals);
        display_mismatches = display_mismatches + nnz(~matches);
        assert(all(matches), ...
            'NODDI reference changed at reported precision in %s.', name);
    end

    observed_p = double(noddi_actual.p_raw(noddi_order));
    expected_p = double(noddi.p_raw);
    observed_p_adj = double(noddi_actual.p_adj(noddi_order));
    expected_p_adj = double(noddi.p_adj);
    assert(all(p_display_matches(observed_p, expected_p)) && ...
           all(p_display_matches(observed_p_adj, expected_p_adj)), ...
        'NODDI p values changed at manuscript precision.');
    adjusted_rows = logical(noddi_actual.Adjusted(noddi_order));
    assert(all(abs(observed_p_adj(adjusted_rows) - ...
        min(observed_p(adjusted_rows) .* 98, 1)) <= absolute_tolerance) && ...
        all(isnan(observed_p_adj(~adjusted_rows))), ...
        'NODDI rows do not contain internally valid m=98 p values.');
    assert(isequal(observed_p_adj < 0.05, expected_p_adj < 0.05), ...
        'NODDI adjusted significance decisions changed.');

    observed_lambda = double(noddi_actual.Lambda(noddi_order));
    expected_lambda = double(noddi.Lambda);
    lambda_difference = abs(observed_lambda - expected_lambda);
    finite_lambda = isfinite(lambda_difference);
    if any(finite_lambda)
        max_lambda_difference = max(lambda_difference(finite_lambda));
        lambda_rows_different = nnz(lambda_difference(finite_lambda) > ...
            absolute_tolerance);
    else
        max_lambda_difference = 0;
        lambda_rows_different = 0;
    end
    p_difference_noddi = abs(observed_p - expected_p);
    p_difference_noddi = p_difference_noddi(isfinite(p_difference_noddi));
    if ~isempty(p_difference_noddi)
        max_noddi_difference = max(max_noddi_difference, ...
            max(p_difference_noddi));
    end

    validation = struct( ...
        'legacy_adjusted_rows_checked', 56, ...
        'legacy_maximum_p_raw_difference', max(p_difference), ...
        'legacy_maximum_p_adj_m98_difference', max(p98_difference), ...
        'legacy_p_display_mismatches', nnz(~legacy_p_display), ...
        'legacy_p_adj_display_mismatches', nnz(~legacy_p98_display), ...
        'noddi_rows_checked', 84, ...
        'noddi_maximum_numeric_difference', max_noddi_difference, ...
        'noddi_display_mismatches', display_mismatches, ...
        'noddi_lambda_rows_different', lambda_rows_different, ...
        'noddi_maximum_lambda_difference', max_lambda_difference, ...
        'absolute_tolerance', absolute_tolerance, ...
        'reference_contract', [ ...
            'keys/cohorts/family exact; AUC/R2/CI 3 decimals; AIC 2 decimals; ' ...
            'p values 4 decimals or <0.0001; significance exact; lambda diagnostic'], ...
        'status', 'PASS');
end


function matches = numeric_display_matches(observed, expected, decimals)
    assert(isequal(size(observed), size(expected)), ...
        'Display-comparison vectors differ in size.');
    matches = false(size(observed));
    both_missing = isnan(observed) & isnan(expected);
    matches(both_missing) = true;
    finite = isfinite(observed) & isfinite(expected);
    format = sprintf('%%.%df', decimals);
    finite_indices = find(finite);
    for ii = 1:numel(finite_indices)
        index = finite_indices(ii);
        matches(index) = strcmp(sprintf(format, observed(index)), ...
            sprintf(format, expected(index)));
    end
end


function matches = p_display_matches(observed, expected)
    assert(isequal(size(observed), size(expected)), ...
        'P-value comparison vectors differ in size.');
    matches = false(size(observed));
    for ii = 1:numel(observed)
        matches(ii) = strcmp(format_p_for_reference(observed(ii)), ...
            format_p_for_reference(expected(ii)));
    end
end


function value = format_p_for_reference(p_value)
    if isnan(p_value)
        value = 'NA';
    elseif p_value < 0.0001
        value = '<0.0001';
    else
        value = sprintf('%.4f', p_value);
    end
end


function values = logical_reference_column(column)
    if islogical(column)
        values = column;
    elseif isnumeric(column)
        assert(all(column == 0 | column == 1), ...
            'Numeric logical reference column must contain only 0/1.');
        values = logical(column);
    else
        text = lower(strtrim(string(column)));
        assert(all(text == "true" | text == "false" | text == "1" | text == "0"), ...
            'Text logical reference column has an invalid value.');
        values = text == "true" | text == "1";
    end
end


function order = match_reference_rows(actual, reference, include_adjusted)
% Map each frozen row to exactly one active row using stable model identity.
    order = zeros(height(reference),1);
    for ri = 1:height(reference)
        match = string(actual.Outcome) == string(reference.Outcome(ri)) & ...
            string(actual.Metric) == string(reference.Metric(ri)) & ...
            string(actual.Model) == string(reference.Model(ri));
        if include_adjusted
            match = match & actual.Adjusted == logical(reference.Adjusted(ri));
        end
        index = find(match);
        assert(isscalar(index), ...
            'Expected one active reference match for %s/%s/%s.', ...
            string(reference.Outcome(ri)), string(reference.Metric(ri)), ...
            string(reference.Model(ri)));
        order(ri) = index;
    end
    assert(numel(unique(order)) == height(reference), ...
        'Reference-to-active row mapping was not one-to-one.');
end


function [full_path, matched_path] = write_scope_figure_model_tables( ...
        full_rows, matched_rows, output_dir, all_metrics_mode)
% Unified mode writes all 70 actively refitted rows. Compatibility mode
% retains the frozen 40 legacy rows and replaces only its 30 NODDI rows.
    if all_metrics_mode
        assert(numel(full_rows) == 70 && numel(matched_rows) == 70, ...
            'Expected 70 full and 70 matched actively refitted rows.');
        full_table = order_figure_table(normalize_figure_table( ...
            struct2table(full_rows)));
        matched_table = order_figure_table(normalize_figure_table( ...
            struct2table(matched_rows)));
        full_path = fullfile(output_dir,'figure6_full_models.csv');
        matched_path = fullfile(output_dir,'figure6_matched_models.csv');
        writetable(full_table,full_path);
        writetable(matched_table,matched_path);
        return;
    end

    assert(numel(full_rows) == 30 && numel(matched_rows) == 30, ...
        'Expected 30 full and 30 matched NODDI continuous rows.');
    matlab_dir = fileparts(mfilename('fullpath'));
    helper_dir = fileparts(matlab_dir);
    package_root = fileparts(helper_dir);
    reference_dir = fullfile(package_root,'results','model_results');
    full_reference = readtable(fullfile(reference_dir,'figure6_full_models.csv'), ...
        'TextType','string');
    matched_reference = readtable(fullfile(reference_dir,'figure6_matched_models.csv'), ...
        'TextType','string');

    full_legacy = keep_legacy_figure_rows(full_reference);
    matched_legacy = keep_legacy_figure_rows(matched_reference);
    assert(height(full_legacy) == 40 && height(matched_legacy) == 40, ...
        'Expected 40 legacy (non-NODDI) rows in each Figure 6 reference table.');

    full_new = normalize_figure_table(struct2table(full_rows));
    matched_new = normalize_figure_table(struct2table(matched_rows));
    full_table = order_figure_table([full_legacy; full_new]);
    matched_table = order_figure_table([matched_legacy; matched_new]);
    assert(height(full_table) == 70 && height(matched_table) == 70, ...
        'Expected 70 rows in each complete Figure 6 model table.');

    full_path = fullfile(output_dir,'figure6_full_models.csv');
    matched_path = fullfile(output_dir,'figure6_matched_models.csv');
    writetable(full_table,full_path);
    writetable(matched_table,matched_path);
end


function T = keep_legacy_figure_rows(T)
    T = normalize_figure_table(T);
    key = string(T.ClinicalMetric);
    is_noddi = endsWith(key,'-NDI') | endsWith(key,'-ODI') | ...
        endsWith(key,'-ISOVF');
    T = T(~is_noddi,:);
end


function T = normalize_figure_table(T)
    expected = {'ClinicalMetric','Model','N','R2','R2_Demo'};
    assert(all(ismember(expected,T.Properties.VariableNames)), ...
        'Figure model table is missing a required column.');
    T = T(:,expected);
    T.ClinicalMetric = string(T.ClinicalMetric);
    T.Model = lower(strtrim(string(T.Model)));
    T.N = double(T.N);
    T.R2 = double(T.R2);
    T.R2_Demo = double(T.R2_Demo);
    assert(all(isfinite(T.N)) && all(isfinite(T.R2)) && all(isfinite(T.R2_Demo)), ...
        'Figure model table contains a non-finite numeric value.');
end


function T = order_figure_table(T)
    clinical = {'T25FW','x9HPTD','x9HPTND','SDMTcorrect','MSFC_SDMT'};
    metrics = {'T1','MTR','FA','MD','NDI','ODI','ISOVF'};
    models = {'tract-based','classic'};
    order = zeros(70,1);
    oi = 0;
    for ci = 1:numel(clinical)
        for mi = 1:numel(metrics)
            key = sprintf('%s-%s',clinical{ci},metrics{mi});
            for fi = 1:numel(models)
                index = find(T.ClinicalMetric == key & T.Model == models{fi});
                assert(isscalar(index), ...
                    'Expected exactly one %s row for %s.', models{fi}, key);
                oi = oi + 1;
                order(oi) = index;
            end
        end
    end
    assert(height(T) == oi, 'Unexpected extra Figure 6 model rows.');
    T = T(order,:);
end


function rows = run_legacy_validation(T_clin, stats_dir, metrics_dir, nrf_preds_base, ...
        ico_regions, bin_measures, bin_thresh, cont_measures, cont_labels, ...
        demo_cols_bin, demo_cols_cont, n_boot, alpha_ci)
% Refit only the FA/MD tract rows used as the manuscript-equivalence gate.
% The case selection and model calls are direct copies of the corresponding
% tract branch in rev_qMRI.m.
    metrics = {'FA','MD'};
    T_ico = readtable(fullfile(stats_dir, 'icometrixLesionStats.xlsx'));
    T_crf_full = outerjoin(T_clin, T_ico, 'Keys','SubjectID','MergeKeys',true);
    T_crf_ms = T_crf_full(~startsWith(T_crf_full.SubjectID, 'sub-C'), :);
    rows = struct([]);

    for mi = 1:numel(metrics)
        metric = metrics{mi};
        generated_path = fullfile(metrics_dir, ...
            sprintf('GroupTract%s_All.xlsx', metric));
        if isfile(generated_path)
            path = generated_path;
        else
            path = fullfile(stats_dir, sprintf('GroupTract%s_All.xlsx', metric));
        end
        require_file(path);
        T_tract = readtable(path);
        T_full = outerjoin(T_clin, T_tract, 'Keys','SubjectID','MergeKeys',true);
        T_nrf_full.(metric) = T_full;
        T_nrf_ms.(metric) = T_full(~startsWith(T_full.SubjectID, 'sub-C'), :);
    end

    lambda_grid_bin = logspace(-6, 6, 50);
    for ci = 1:numel(bin_measures)
        measure = bin_measures{ci};
        threshold = bin_thresh(ci);
        for mi = 1:numel(metrics)
            metric = metrics{mi};
            nrf_preds = validate_preds(nrf_preds_base, T_nrf_ms.(metric), ...
                sprintf('NRF-%s', metric));
            crf_preds = validate_preds(build_crf_preds(ico_regions, metric), ...
                T_crf_ms, sprintf('CRF-%s', metric));
            for adjusted = [false true]
                if adjusted, demo_cols = demo_cols_bin; else, demo_cols = {}; end
                [X, ~, y, Demo, n] = get_joint_cases( ...
                    T_nrf_ms.(metric), T_crf_ms, nrf_preds, crf_preds, ...
                    measure, demo_cols, true, threshold, false, true, 'binary');
                rng(42);
                [LP, lambda] = fit_ridge_lp_bin(zscore(X), y, ...
                    lambda_grid_bin, sprintf('%s|Tract-based %s|validation', ...
                    measure, metric));
                res = eval_binary(LP, y, Demo, n, n_boot, alpha_ci, ...
                    sprintf('%s|Tract-based %s|validation', measure, metric));
                rows = append_canonical_row(rows, measure, measure, metric, ...
                    'Tract-based', adjusted, demo_cols, ...
                    'legacy joint tract+classical completeness (validation only)', ...
                    res, lambda, true, y, NaN);
            end
        end
    end

    for ci = 1:numel(cont_measures)
        measure = cont_measures{ci};
        label = cont_labels{ci};
        for mi = 1:numel(metrics)
            metric = metrics{mi};
            nrf_preds = validate_preds(nrf_preds_base, T_nrf_full.(metric), ...
                sprintf('NRF-%s', metric));
            for adjusted = [false true]
                if adjusted, demo_cols = demo_cols_cont; else, demo_cols = {}; end
                [X, y, Demo, n] = get_complete_cases( ...
                    T_nrf_full.(metric), nrf_preds, measure, demo_cols, ...
                    false, [], [], true, 'log');
                rng(42);
                res = run_model_cont(X, y, Demo, n, ...
                    sprintf('%s|Tract-based %s|validation', label, metric));
                rows = append_canonical_row(rows, label, measure, metric, ...
                    'Tract-based', adjusted, demo_cols, ...
                    'independent tract complete cases', ...
                    res, res.lambda, false, [], NaN);
            end
        end
    end
end


function text = family_definition(family_size)
assert(family_size == 98, ...
    'This frozen reviewer package requires globalFamilySize = 98.');
text = ['98 = existing 56 +Demo qMRI tests plus 42 +Demo NODDI tests ' ...
    '(3 metrics x 7 outcomes x 2 model frameworks)'];
end
function values = append_struct(values, item)
if isempty(values), values = item; else, values(end+1) = item; end
end


function require_file(path)
assert(isfile(path), 'Required file not found: %s', path);
end


function filename = get_filename(path)
[~,stem,extension] = fileparts(path);
filename = [stem extension];
end
