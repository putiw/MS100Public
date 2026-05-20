function results = rev_corr(metric_type, varargin)
% Optional output: struct array, one element per heatmap produced.
%   results(i).corr_matrix     — absolute correlation values
%   results(i).p_raw           — raw (uncorrected) p-values
%   results(i).display_names   — row labels
%   results(i).clinical_labels — column labels
%   results(i).title           — figure title string
% REV_CORR  Updated correlation script for Brain Communications revision.
%
% Addresses Reviewer 1, Major Point 2: T25FW, 9HPTD, and 9HPTND are
% log-transformed before computing Pearson correlations, consistent with
% the regression models in paper_multi_continuous_ridge.m. EDSS is left
% on its original scale (not log-transformed) as it is ordinal.
%
% Usage (mirrors existing f10_corr_singleTract.m interface):
%   rev_corr('volume')                        % Lesion volume
%   rev_corr('norm')                          % Normalized lesion volume
%   rev_corr('number')                        % Number of lesions
%   rev_corr('T1')                            % T1 (all tissue)
%   rev_corr('FA', 'tissue', 'Lesion')        % FA in lesions
%   rev_corr('MD', 'tissue', 'NAWM')          % MD in NAWM
%   rev_corr('MTR', 'tissue', 'All')          % MTR all tissue
%
% Optional parameters:
%   'tissue'    - 'All' (default), 'Lesion', or 'NAWM'
%   'tractType' - 'dwi' (default) or 'template'
%   'groupOnly' - true to plot only grouped tract clusters (not individual
%                 tracts), false (default) to plot individual tracts also

    %% Parse inputs
    % Detect optional trailing title string, e.g. rev_corr(..., 'My Title')
    custom_title = '';
    known_keys = {'tissue', 'tractType', 'groupOnly', 'plot'};
    if ~isempty(varargin) && ischar(varargin{end}) && ...
            ~ismember(varargin{end}, known_keys) && ...
            ~ismember(varargin{end}, {'All','Lesion','NAWM','dwi','template'}) && ...
            mod(numel(varargin), 2) == 1
        custom_title = varargin{end};
        varargin     = varargin(1:end-1);
    end

    p = inputParser;
    addRequired(p, 'metric_type', @ischar);
    addParameter(p, 'tissue',    'All',  @(x) ismember(x, {'All','Lesion','NAWM','Tail'}));
    addParameter(p, 'tractType', 'dwi',  @(x) ismember(x, {'dwi','template'}));
    addParameter(p, 'groupOnly', false,  @islogical);
    addParameter(p, 'plot',      true,   @islogical);
    parse(p, metric_type, varargin{:});

    metric_type = p.Results.metric_type;
    tissue_type = p.Results.tissue;
    tract_type  = p.Results.tractType;
    group_only  = p.Results.groupOnly;
    do_plot     = p.Results.plot;

    %% Initialise output
    results = struct('corr_matrix',{}, 'p_raw',{}, 'display_names',{}, ...
                     'clinical_labels',{}, 'title',{});

    %% Setup
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));
    cfg = load_config();

    % Clinical scores included in the correlation analysis.
    % T25FW, x9HPTD, x9HPTND will be log-transformed (see apply_log_transform).
    % EDSS is kept on the original ordinal scale.
    clinical_scores = {'EDSS', 'T25FW', 'x9HPTD', 'x9HPTND' , 'MSPro'};

    % Display labels for the x-axis of heatmaps (after transformation)
    clinical_labels = {'EDSS', 'log(T25FW)', 'log(9HPT-D)', 'log(9HPT-ND)', 'MSPro'};

    %% Handle WB (whole brain) special case — early return
    % Usage: rev_corr('WB') or rev_corr('WB', 'My Title')
    % Shows a single heatmap with WBLN and WBLV vs all clinical outcomes.
    % WBLnorm is not stored in GroupTractLesionLoad.xlsx and is therefore omitted.
    if strcmp(metric_type, 'WB')
        baseDir = fullfile(cfg.bidsDir, cfg.statsDir);
        data_wb = readtable(fullfile(baseDir, 'clinicalScore.xlsx'));
        Ttract  = readtable(fullfile(baseDir, 'GroupTractLesionLoad.xlsx'));
        data_wb = outerjoin(data_wb, Ttract, 'Keys', 'SubjectID', 'MergeKeys', true);

        % Same clinical preprocessing as all other metric types
        data_wb = clean_clinical_data(data_wb, clinical_scores);
        for fld = {'x9HPTD', 'x9HPTND'}
            if ismember(fld{1}, data_wb.Properties.VariableNames)
                col = data_wb.(fld{1});
                col(col <= 0) = NaN;
                data_wb.(fld{1}) = log(col);
            end
        end
        data_wb = remove_outliers(data_wb, clinical_scores);

        % Find available WB columns (WBLN, WBLV; WBLnorm not present)
        all_vars_wb    = data_wb.Properties.VariableNames;
        wb_candidates  = {'WBLN',              'WBLV'};
        wb_labels_all  = {'WB Lesion Number',  'WB Lesion Volume'};
        wb_metrics = {};
        wb_labels  = {};
        for wi = 1:numel(wb_candidates)
            if ismember(wb_candidates{wi}, all_vars_wb)
                wb_metrics{end+1} = wb_candidates{wi}; %#ok<AGROW>
                wb_labels{end+1}  = wb_labels_all{wi}; %#ok<AGROW>
            end
        end

        if isempty(wb_metrics)
            warning('rev_corr:WB', 'No WB columns (WBLN/WBLV) found in GroupTractLesionLoad.xlsx');
            return;
        end

        title_str = 'Whole Brain Lesion Load';
        if ~isempty(custom_title), title_str = custom_title; end

        [cm, pr, dn] = plot_correlations(data_wb, wb_metrics, clinical_scores, ...
                          clinical_labels, title_str, '', wb_labels, do_plot);
        results(end+1) = make_entry(cm, pr, dn, clinical_labels, title_str); %#ok<AGROW>
        return;
    end

    %% Determine metric class
    lesion_load_metrics  = {'norm', 'volume', 'number'};
    quantitative_metrics = {'T1', 'FA', 'MD', 'MTR', 'AD', 'RD'};

    if ismember(metric_type, lesion_load_metrics)
        is_lesion_load = true;
    elseif ismember(metric_type, quantitative_metrics)
        is_lesion_load = false;
    else
        error('Invalid metric_type. Choose from: WB, %s, or %s', ...
            strjoin(lesion_load_metrics, ', '), strjoin(quantitative_metrics, ', '));
    end

    %% Load Data — use IDENTICAL pipeline as arab_health_figure.m
    if is_lesion_load
        % For lesion load metrics, load grouped data directly
        baseDir = fullfile(cfg.bidsDir, cfg.statsDir);
        data = readtable(fullfile(baseDir, 'clinicalScore.xlsx'));
        if group_only
            Ttract = readtable(fullfile(baseDir, 'GroupTractLesionLoad.xlsx'));
        else
            Ttract = readtable(fullfile(baseDir, 'TractLesionLoad.xlsx'));
        end
        data = outerjoin(data, Ttract, 'Keys', 'SubjectID', 'MergeKeys', true);
        suffix       = get_lesion_suffix(metric_type);
        title_prefix = ['Lesion ' metric_type];
    else
        % For quantitative metrics, prefer a tissue-specific file
        % (e.g. GroupTractMTR_All.xlsx) and fall back to the combined file
        % (e.g. GroupTractT1.xlsx) whose columns carry _Lesion/_NAWM suffixes.
        baseDir = fullfile(cfg.bidsDir, cfg.statsDir);
        Tclin   = readtable(fullfile(baseDir, 'clinicalScore.xlsx'));
        Tico    = readtable(fullfile(baseDir, 'icometrixLesionStats.xlsx'));
        data    = outerjoin(Tclin, Tico, 'Keys', 'SubjectID', 'MergeKeys', true);
        metric_file = fullfile(baseDir, sprintf('GroupTract%s_All.xlsx', metric_type));
        Tmetric = readtable(metric_file);
        data    = outerjoin(data, Tmetric, 'Keys', 'SubjectID', 'MergeKeys', true);
        suffix       = ['_' tissue_type];
        title_prefix = [metric_type ' (' tissue_type ')'];
    end

    %% Preprocess Clinical Data — IDENTICAL to arab_health_figure.m
    % 1. clean_clinical_data: log-transforms T25FW & winsorizes (same helper function)
    data = clean_clinical_data(data, clinical_scores);
    
    % 2. ADDITIONALLY log-transform 9HPTD and 9HPTND (Reviewer 1 Major Point 2)
    %    clean_clinical_data only handles T25FW; the reviewer asked for all timed metrics
    for fld = {'x9HPTD', 'x9HPTND'}
        if ismember(fld{1}, data.Properties.VariableNames)
            col = data.(fld{1});
            col(col <= 0) = NaN;
            data.(fld{1}) = log(col);
        end
    end
    
    % 3. Remove outliers (z > 3), same as arab_health_figure
    data = remove_outliers(data, clinical_scores);

    % 4. For _Lesion columns, treat 0 as NaN: a value of 0 means the subject
    %    has no lesion voxels in that tract (e.g. controls), not a real measurement.
    %    This does not affect _Tail, _All, or _NAWM columns.
    if strcmp(tissue_type, 'Lesion')
        lesion_cols = data.Properties.VariableNames(endsWith(data.Properties.VariableNames, '_Lesion'));
        for ci = 1:numel(lesion_cols)
            col = data.(lesion_cols{ci});
            col(col == 0) = NaN;
            data.(lesion_cols{ci}) = col;
        end
    end

    %% Get tract / metric columns
    all_vars = data.Properties.VariableNames;
    tract_metrics = all_vars(endsWith(all_vars, suffix));
    tract_metrics = setdiff(tract_metrics, [{'SubjectID'}, clinical_scores]);

    % When showing grouped data, remove the L/R lateralized columns
    % and keep only the combined (bilateral) ones
    if group_only && ~is_lesion_load
        tract_metrics = tract_metrics(~contains(tract_metrics, 'L_') & ...
                                       ~contains(tract_metrics, 'R_'));
    elseif group_only && is_lesion_load
        % For lesion load, L/R are encoded as e.g. AssociationLLV, AssociationRLV
        % Keep only ones without L/R before the suffix
        keep = true(size(tract_metrics));
        for i = 1:numel(tract_metrics)
            name = tract_metrics{i};
            before_suffix = strrep(name, suffix, '');
            % If the name ends with L or R before suffix, it's lateralized
            if length(before_suffix) > 1 && (before_suffix(end) == 'L' || before_suffix(end) == 'R')
                keep(i) = false;
            end
        end
        tract_metrics = tract_metrics(keep);
    end

    % Remove WB (whole white matter) row from tract-based table
    tract_metrics = tract_metrics(~startsWith(tract_metrics, 'WB'));

    %% Get classical (icometrix) metric columns
    if ~is_lesion_load
        % For quantitative metrics: periventricularT1, juxtacorticalT1, etc.
        ico_names = {'periventricular', 'juxtacortical', 'infratentorial', 'deepwhitematter'};
        ico_metrics = {};
        for i = 1:numel(ico_names)
            col_name = [ico_names{i} metric_type];
            if ismember(col_name, all_vars)
                ico_metrics{end+1} = col_name;
            end
        end
        ico_suffix = metric_type;  % suffix for display name stripping
    else
        % For lesion load: periventricularLV, juxtacorticalLnorm, etc.
        ico_names = {'periventricular', 'juxtacortical', 'infratentorial', 'deepwhitematter'};
        % Load icometrix lesion load if not already in the table
        baseDir2 = fullfile(cfg.bidsDir, cfg.statsDir);
        if ~ismember(['periventricular' suffix], all_vars)
            Tico = readtable(fullfile(baseDir2, 'icometrixLesionLoad.xlsx'));
            data = outerjoin(data, Tico, 'Keys', 'SubjectID', 'MergeKeys', true);
            all_vars = data.Properties.VariableNames;
        end
        ico_metrics = {};
        for i = 1:numel(ico_names)
            col_name = [ico_names{i} suffix];
            if ismember(col_name, all_vars)
                ico_metrics{end+1} = col_name;
            end
        end
        ico_suffix = suffix;
    end

    %% Plot Figure 1: Tract-based
    if ~isempty(tract_metrics)
        base = custom_title; if isempty(base), base = title_prefix; end
        t1 = ['Tract-Based ' base];
        [cm, pr, dn] = plot_correlations(data, tract_metrics, clinical_scores, ...
                          clinical_labels, t1, suffix, [], do_plot);
        results(end+1) = make_entry(cm, pr, dn, clinical_labels, t1); %#ok<AGROW>
    else
        warning('No tract metrics found for the specified criteria.');
    end

    %% Plot Figure 2: Classical (icometrix)
    if ~isempty(ico_metrics) && ~strcmp(tissue_type, 'NAWM')
        base = custom_title; if isempty(base), base = title_prefix; end
        t2 = ['Classical ' base];
        [cm, pr, dn] = plot_correlations(data, ico_metrics, clinical_scores, ...
                          clinical_labels, t2, ico_suffix, [], do_plot);
        results(end+1) = make_entry(cm, pr, dn, clinical_labels, t2); %#ok<AGROW>
    end
end


%% =========================================================================
%  Data loading
%  =========================================================================

function data = load_lesion_load_data(cfg, clinical_scores, metric_type, group_only)
    baseDir = fullfile(cfg.bidsDir, cfg.statsDir);
    Tclin   = readtable(fullfile(baseDir, 'clinicalScore.xlsx'));
    if group_only
        Ttract  = readtable(fullfile(baseDir, 'GroupTractLesionLoad.xlsx'));
    else
        Ttract  = readtable(fullfile(baseDir, 'TractLesionLoad.xlsx'));
    end
    data    = outerjoin(Tclin, Ttract, 'Keys', 'SubjectID', 'MergeKeys', true);
    data    = clean_missing(data, clinical_scores);
end

function data = load_quantitative_data(cfg, clinical_scores, metric_type, tract_type, group_only)
    baseDir = fullfile(cfg.bidsDir, cfg.statsDir);
    Tclin   = readtable(fullfile(baseDir, 'clinicalScore.xlsx'));

    if group_only
        filename = ['GroupTract', upper(metric_type), '.xlsx'];
    else
        if strcmp(tract_type, 'template')
            filename = cfg.excelOutputs.template.individualStats.(lower(metric_type));
        else
            filename = cfg.excelOutputs.dwi.individualStats.(lower(metric_type));
        end
    end

    metric_file = fullfile(baseDir, filename);
    if ~exist(metric_file, 'file')
        error('Metric file not found: %s', metric_file);
    end
    Tmetric = readtable(metric_file);
    data    = outerjoin(Tclin, Tmetric, 'Keys', 'SubjectID', 'MergeKeys', true);
    data    = clean_missing(data, clinical_scores);
end

function data = clean_missing(data, clinical_scores)
    % Remove control subjects (no clinical scores expected)
    if ismember('SubjectID', data.Properties.VariableNames)
        control_mask = startsWith(data.SubjectID, 'sub-C');
        data = data(~control_mask, :);
    end
end

function suffix = get_lesion_suffix(metric_type)
    switch metric_type
        case 'norm',   suffix = 'Lnorm';
        case 'volume', suffix = 'LV';
        case 'number', suffix = 'LN';
    end
end


%% =========================================================================
%  Log transformation
%  =========================================================================

function data = apply_log_transform(data)
% APPLY_LOG_TRANSFORM  Log-transforms timed motor scores before correlation.
%
% Rationale (Reviewer 1, Major Point 2):
%   T25FW, 9HPT-D, and 9HPT-ND are positively skewed timed measures.
%   Natural log transformation (consistent with paper_multi_continuous_ridge.m)
%   produces approximately normal distributions suitable for Pearson correlation.
%   EDSS is ordinal and is NOT transformed here.
%
%   Subjects with zero or missing values are set to NaN to exclude them
%   from the pairwise correlation rather than producing -Inf on log(0).

    timed_scores = {'T25FW', 'x9HPTD', 'x9HPTND'};

    for i = 1:numel(timed_scores)
        fld = timed_scores{i};
        if ismember(fld, data.Properties.VariableNames)
            col = data.(fld);
            % Replace non-positive with NaN before log (avoids log(0) = -Inf)
            col(col <= 0) = NaN;
            data.(fld) = log(col);
            fprintf('Log-transformed %s: range [%.3f, %.3f], n=%d valid\n', ...
                fld, min(data.(fld), [], 'omitnan'), ...
                max(data.(fld), [], 'omitnan'), sum(~isnan(data.(fld))));
        end
    end
end


%% =========================================================================
%  Correlation computation
%  =========================================================================

function [corr_matrix, p_matrix] = compute_correlations(data, metrics, clinical_scores)
% COMPUTE_CORRELATIONS  Pairwise partial correlations with pairwise-complete obs.
%
%   Zeros in clinical scores are treated as missing.
%   Uses Spearman for ordinal variables (EDSS, MSPro) and Pearson for continuous.
%   Controls for Age and Gender via partial correlation.

    num_metrics = numel(metrics);
    num_scores  = numel(clinical_scores);
    corr_matrix = NaN(num_metrics, num_scores);
    p_matrix    = NaN(num_metrics, num_scores);
    
    % Covariates for partial correlation.
    % When healthy controls are included (continuous outcomes), disease
    % duration is omitted because controls lack this variable.
    % When only MS patients are analysed (binary outcomes such as EDSS,
    % MSPro), disease duration is added as a third covariate.
    % For correlations the cohort is MS-only, so all three are used.
    confounds = {'Age', 'Gender', 'DurationOfDisease'};

    for k = 1:num_metrics
        for j = 1:num_scores
            x = data.(metrics{k});
            y = data.(clinical_scores{j});
            score_name = clinical_scores{j};

            % Pairwise complete: exclude NaN
            valid = ~isnan(x) & ~isnan(y);
            
            % Only exclude zeros for timed measures (not EDSS)
            if ismember(score_name, {'T25FW', 'x9HPTD', 'x9HPTND'})
                valid = valid & y ~= 0;
            end
            
            % Also require valid covariates
            for c = 1:length(confounds)
                if ismember(confounds{c}, data.Properties.VariableNames)
                    c_data = data.(confounds{c});
                    if iscategorical(c_data) || isstring(c_data) || iscell(c_data)
                        valid = valid & ~isundefined(categorical(c_data));
                    else
                        valid = valid & ~isnan(c_data);
                    end
                end
            end

            if sum(valid) < 3
                continue;   % Leave as NaN
            end
            
            x_valid = x(valid);
            y_valid = y(valid);
            
            % Build confound matrix
            confound_matrix = [];
            for c = 1:length(confounds)
                if ismember(confounds{c}, data.Properties.VariableNames)
                    c_data = data.(confounds{c})(valid);
                    if iscategorical(c_data) || isstring(c_data) || iscell(c_data)
                        [~, ~, c_numeric] = unique(cellstr(c_data));
                        confound_matrix = [confound_matrix, c_numeric(:)];
                    else
                        confound_matrix = [confound_matrix, c_data(:)];
                    end
                end
            end

            if ismember(score_name, {'EDSS', 'MSPro'})
                % Partial Spearman: MATLAB built-in ranks all variables
                % including confounds before computing partial correlation
                if isempty(confound_matrix)
                    [rmat, pmat] = corr([x_valid, y_valid], 'Type', 'Spearman');
                else
                    [rmat, pmat] = partialcorr([x_valid, y_valid], confound_matrix, ...
                                               'Type', 'Spearman');
                end
                r    = rmat(1, 2);
                pval = pmat(1, 2);
            else
                % Partial Pearson: MATLAB built-in
                if isempty(confound_matrix)
                    [rmat, pmat] = corr([x_valid, y_valid], 'Type', 'Pearson');
                else
                    [rmat, pmat] = partialcorr([x_valid, y_valid], confound_matrix, ...
                                               'Type', 'Pearson');
                end
                r    = rmat(1, 2);
                pval = pmat(1, 2);
            end
            
            corr_matrix(k, j) = abs(r);
            p_matrix(k, j)    = pval;
        end
    end
end


%% =========================================================================
%  Plotting
%  =========================================================================

function [corr_matrix, p_raw, display_names] = plot_correlations(data, tract_metrics, ...
                            clinical_scores, clinical_labels, title_prefix, suffix, ...
                            display_names_in, do_plot)
% Returns raw (uncorrected) corr_matrix, p_raw, and display_names.
% do_plot (default true): set false to suppress figure generation.
if nargin < 8, do_plot = true; end

    % Build display names by stripping suffix (or use caller-supplied names)
    if nargin >= 7 && ~isempty(display_names_in)
        display_names = display_names_in(:);
    else
    display_names = cell(size(tract_metrics));
    for i = 1:numel(tract_metrics)
        name = tract_metrics{i};
        
        % For grouped lesion load files, suffixes lack the underscore and have L/R
        % e.g. AssociationRLV or AssociationLN
        name = strrep(name, suffix, '');
        name = strrep(name, strrep(suffix, '_', ''), ''); 
        name = strrep(name, 'LLN', ' (L)'); name = strrep(name, 'RLN', ' (R)'); name = strrep(name, 'LN', '');
        name = strrep(name, 'LLV', ' (L)'); name = strrep(name, 'RLV', ' (R)'); name = strrep(name, 'LV', '');
        name = strrep(name, 'LLnorm', ' (L)'); name = strrep(name, 'RLnorm', ' (R)'); name = strrep(name, 'Lnorm', '');
        
        name = strrep(name, '_', ' ');
        name = regexprep(name, '([a-z])([A-Z])', '$1 $2');

        % Expand known abbreviations / fix casing
        name = strrep(name, 'ProjectionBrainstem', 'Projection Brainstem');
        if strcmp(strtrim(name), 'PB'), name = 'Projection Brainstem'; end

        % Fix icometrix region label casing
        name = strrep(name, 'deepwhitematter', 'Deep White Matter');
        name = strrep(name, 'periventricular',  'Periventricular');
        name = strrep(name, 'juxtacortical',    'Juxtacortical');
        name = strrep(name, 'infratentorial',   'Infratentorial');

        display_names{i} = strtrim(name);
    end
    end  % end else (display_names_in was empty)

    [corr_matrix, p_matrix] = compute_correlations(data, tract_metrics, clinical_scores);
    p_raw = p_matrix;   % save before correction

    % Note: when called from result252.m, correction is applied externally
    % (Bonferroni across all 130 tests). Raw p-values are returned in
    % results.p_raw for that purpose.

    if ~do_plot, return; end

    % Single-column paper figure (two-column layout, ~3.5 inch column width)
    num_tracts = numel(tract_metrics);
    num_scores = numel(clinical_scores);
    fig_w_in  = 3.5;                              % fixed single-column width (inches)
    cell_h_in = 0.20;                             % height per tract row (inches)
    margin_h_in = 1.1;                            % title + x-labels + colorbar
    fig_h_in  = margin_h_in + num_tracts * cell_h_in;
    figure('Units', 'inches', 'Position', [1, 1, fig_w_in, fig_h_in], ...
           'PaperUnits', 'inches', 'PaperSize', [fig_w_in, fig_h_in], ...
           'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial');

    imagesc(corr_matrix);
    colorbar;
    colormap(flip(hot));
    caxis([0 1]);

    title([title_prefix ' Correlations'], ...
          'FontSize', 7, 'FontWeight', 'bold');

    ax_fs = 6;
    set(gca, 'XTick', 1:num_scores, 'XTickLabel', clinical_labels, 'FontSize', ax_fs);
    set(gca, 'YTick', 1:num_tracts, 'YTickLabel', display_names,   'FontSize', ax_fs);
    xtickangle(30);

    % Annotate cells: bold + * for significant
    cell_fs = 6;
    for k = 1:num_tracts
        for j = 1:num_scores
            if ~isnan(corr_matrix(k, j))
                sig = p_matrix(k, j) < 0.05;
                txt = sprintf('%.2f', corr_matrix(k, j));
                if sig, txt = [txt '*']; end %#ok<AGROW>
                if sig, fw = 'bold'; else, fw = 'normal'; end
                if corr_matrix(k, j) > 0.6, fc = 'w'; else, fc = 'k'; end
                text(j, k, txt, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment',   'middle', ...
                    'FontSize',   cell_fs, ...
                    'FontWeight', fw, ...
                    'Color', fc);
            end
        end
    end
end


%% =========================================================================
%  Data Cleaning Helpers
%  =========================================================================

function s = make_entry(cm, pr, dn, cl, title_str)
% Pack one heatmap's raw data into a struct for external joint correction.
    s = struct('corr_matrix', cm, 'p_raw', pr, ...
               'display_names', {dn}, 'clinical_labels', {cl}, 'title', title_str);
end


function data = apply_winsorization(data, clinical_scores)
    for ii = 1:numel(clinical_scores)
        fld = clinical_scores{ii};
        if ismember(fld, data.Properties.VariableNames)
            col = data.(fld);
            m   = mean(col,'omitnan');
            sd  = std(col,'omitnan');
            hi  = m + 2.5*sd;
            lo  = m - 2.5*sd;
            % clamp extreme values
            data.(fld)(col > hi) = hi;
            data.(fld)(col < lo) = lo;
        end
    end
end

function data = remove_outliers(data, clinical_scores)
    threshold = 3;
    for i = 1:length(clinical_scores)
        score = clinical_scores{i};
        if ~ismember(score, data.Properties.VariableNames)
            continue;
        end
        
        values = data.(score);
        if ~isnumeric(values)
            continue;
        end
        
        valid_idx = ~isnan(values);
        z_scores = abs((values - mean(values, 'omitnan')) / std(values, 'omitnan'));
        outlier_idx = z_scores > threshold;
        
        n_outliers = sum(outlier_idx & valid_idx);
        if n_outliers > 0
            fprintf('  %s: %d outliers removed (z > %g)\n', score, n_outliers, threshold);
            data.(score)(outlier_idx) = NaN;
        end
    end
end

