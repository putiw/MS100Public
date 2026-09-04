function [files, qc, validation] = generate_all_tract_metric_tables(varargin)
% GENERATE_ALL_TRACT_METRIC_TABLES Build every manuscript tract-metric table.
%
% This function is a deliberately small extension of paper_gen_excel_stats.m.
% It preserves that script's four tract groups, L/R/combined masks, lesion /
% NAWM / all-tissue means, MATLAB prctile implementation, and legacy handling
% of controls and missing ROIs for T1, MTR, FA, and MD.  It adds NDI, ODI,
% and AMICO FWF (reported in the manuscript as ISOVF).
%
% NODDI-specific rules:
%   * only voxels inside the transformed AMICO fitting support are sampled;
%   * a fitted value of exactly zero is retained;
%   * voxel-map QC permits at most a 1e-6 numerical bound excursion, with no
%     clipping; every reported aggregate must lie strictly in [0,1];
%   * an empty ROI is NaN, never a fabricated zero;
%   * tails are NDI P10, ODI P90, and FWF/ISOVF P90.
%
% Required name-value inputs:
%   SourceRoot    - BIDS root (or read-only snapshot) containing derivatives/
%   NODDIMapsDir  - transformed_maps root containing sub-* T1-space maps
%   OutputDir     - fresh local directory for aggregate tables and QC
%
% Optional name-value inputs:
%   Metrics            - default: T1,MTR,FA,MD,NDI,ODI,ISOVF
%   Subjects           - explicit subject subset (primarily for validation)
%   LegacyStatsDir     - default: SourceRoot/derivatives/derivativesStats
%   ReferenceNODDIDir  - optional Python-table directory for migration check
%   GeometryTolerance  - maximum affine element difference (default 1e-4)
%   ExpectedSourceSubjects - full-run count gate (default 132; 0 disables)
%   ExpectedNODDISubjects  - full-run count gate (default 132; 0 disables)
%   Force              - overwrite existing generated files (default false)
%   WriteQC            - write two compact validation CSVs (default true)
%
% Outputs retain the filenames expected by the manuscript statistics code.
% Legacy metrics are written as both XLSX and CSV; NODDI metrics as CSV.
% Model-facing CSVs use 17 significant digits and are read back before the
% function returns, so serialization cannot silently perturb a MATLAB double.

p = inputParser;
addParameter(p, 'SourceRoot', '', @is_text_scalar);
addParameter(p, 'NODDIMapsDir', '', @is_text_scalar);
addParameter(p, 'OutputDir', '', @is_text_scalar);
addParameter(p, 'Metrics', {'T1','MTR','FA','MD','NDI','ODI','ISOVF'}, ...
    @(x) iscellstr(x) || isstring(x));
addParameter(p, 'Subjects', {}, @(x) iscellstr(x) || isstring(x));
addParameter(p, 'LegacyStatsDir', '', @is_text_scalar);
addParameter(p, 'ReferenceNODDIDir', '', @is_text_scalar);
addParameter(p, 'GeometryTolerance', 1e-4, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'ValueTolerance', 1e-6, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 0 && x <= 1e-6);
addParameter(p, 'ValidationAbsoluteTolerance', 1e-10, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ValidationRelativeTolerance', 1e-6, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ExpectedSourceSubjects', 132, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ExpectedNODDISubjects', 132, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'Force', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'WriteQC', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

source_root = require_dir(char(p.Results.SourceRoot), 'SourceRoot');
output_dir = char(p.Results.OutputDir);
if isempty(output_dir)
    error('generate_all_tract_metric_tables:MissingOutputDir', ...
        'OutputDir is required.');
end
noddi_maps_dir = char(p.Results.NODDIMapsDir);
metrics = canonical_metrics(p.Results.Metrics);
has_noddi = any(ismember(metrics, {'NDI','ODI','FWF'}));
if has_noddi
    noddi_maps_dir = require_dir(noddi_maps_dir, 'NODDIMapsDir');
end

legacy_stats_dir = char(p.Results.LegacyStatsDir);
if isempty(legacy_stats_dir)
    legacy_stats_dir = fullfile(source_root, 'derivatives', 'derivativesStats');
end
if ~isempty(legacy_stats_dir) && ~isfolder(legacy_stats_dir)
    legacy_stats_dir = '';
end
reference_noddi_dir = char(p.Results.ReferenceNODDIDir);
if ~isempty(reference_noddi_dir)
    reference_noddi_dir = require_dir(reference_noddi_dir, 'ReferenceNODDIDir');
end

explicit_subjects = ~isempty(p.Results.Subjects);
if explicit_subjects
    subject_ids = normalize_subjects(p.Results.Subjects);
else
    subject_ids = discover_subjects(source_root, legacy_stats_dir);
end
if isempty(subject_ids)
    error('generate_all_tract_metric_tables:NoSubjects', ...
        'No sub-* directories were found in the source tree.');
end
if ~explicit_subjects && p.Results.ExpectedSourceSubjects > 0 && ...
        numel(subject_ids) ~= p.Results.ExpectedSourceSubjects
    error('generate_all_tract_metric_tables:UnexpectedSourceCohort', ...
        'Expected %d source subjects; observed %d.', ...
        p.Results.ExpectedSourceSubjects, numel(subject_ids));
end

if ~isfolder(output_dir)
    mkdir(output_dir);
end
planned_outputs = {};
for mi = 1:numel(metrics)
    planned_outputs{end+1} = fullfile(output_dir, ... %#ok<AGROW>
        ['GroupTract' metrics{mi} '_All.csv']);
    if ismember(metrics{mi}, {'T1','MTR','FA','MD'})
        planned_outputs{end+1} = fullfile(output_dir, ... %#ok<AGROW>
            ['GroupTract' metrics{mi} '_All.xlsx']);
    end
end
if p.Results.WriteQC
    planned_outputs = [planned_outputs, { ...
        fullfile(output_dir, 'metric_table_generation_QC.csv'), ...
        fullfile(output_dir, 'metric_table_reference_validation.csv')}];
end
if ~p.Results.Force
    existing = planned_outputs(cellfun(@isfile, planned_outputs));
    if ~isempty(existing)
        error('generate_all_tract_metric_tables:OutputExists', ...
            'Refusing to replace existing output (use Force=true): %s', existing{1});
    end
end
columns = tract_columns();
metric_rows = struct();
for mi = 1:numel(metrics)
    metric_rows.(metrics{mi}) = cell(0, 1);
end
qc_records = empty_qc_records();

fprintf('Generating tract tables for %d subjects and %d metrics.\n', ...
    numel(subject_ids), numel(metrics));

for si = 1:numel(subject_ids)
    subject = subject_ids{si};
    reference_path = find_reference_path(source_root, subject);
    reference_info = read_header_checked(reference_path, ...
        [subject ' T1-grid reference'], p.Results.GeometryTolerance);

    [masks, mask_summary] = load_legacy_masks(source_root, subject, ...
        reference_info, p.Results.GeometryTolerance);

    need_strict_lesion = has_noddi && mask_summary.Complete;
    [lesion_mask, lesion_status] = load_lesion_mask(source_root, subject, ...
        reference_info, need_strict_lesion, p.Results.GeometryTolerance);

    support_mask = [];
    support_info = [];
    support_voxels = NaN;
    if need_strict_lesion
        support_path = fullfile(noddi_maps_dir, subject, ...
            [subject '_space-individual_NODDI-support.nii.gz']);
        [support_data, support_info] = read_volume_checked(support_path, ...
            [subject ' NODDI support'], reference_info, ...
            p.Results.GeometryTolerance);
        if any(~isfinite(double(support_data(:))))
            error('generate_all_tract_metric_tables:InvalidSupport', ...
                '%s support contains non-finite values.', subject);
        end
        support_mask = support_data > 0;
        support_voxels = nnz(support_mask);
        if support_voxels == 0
            error('generate_all_tract_metric_tables:EmptySupport', ...
                '%s has an empty transformed AMICO fitting support.', subject);
        end
    end

    for mi = 1:numel(metrics)
        metric = metrics{mi};
        is_noddi = ismember(metric, {'NDI','ODI','FWF'});

        if is_noddi && ~mask_summary.Complete
            qc_records(end+1) = make_qc_record(subject, metric, ... %#ok<AGROW>
                'SKIPPED_INCOMPLETE_TRACT_MASKS', 0, mask_summary, ...
                NaN, NaN, lesion_status, NaN);
            continue;
        end

        map_path = metric_map_path(source_root, noddi_maps_dir, subject, metric);
        if ~isfile(map_path)
            if is_noddi
                error('generate_all_tract_metric_tables:MissingNODDIMap', ...
                    '%s %s map not found: %s', subject, metric, map_path);
            end
            warning('generate_all_tract_metric_tables:MissingLegacyMap', ...
                '%s %s map not found; preserving legacy missing-row behavior.', ...
                subject, metric);
            values = nan(1, numel(columns));
            metric_rows.(metric){end+1,1} = make_row(subject, columns, values);
            qc_records(end+1) = make_qc_record(subject, metric, ... %#ok<AGROW>
                'MISSING_MAP_LEGACY_NAN_ROW', 1, mask_summary, ...
                NaN, NaN, lesion_status, NaN);
            continue;
        end

        [metric_map, metric_info, map_affine_error] = read_volume_checked( ...
            map_path, [subject ' ' metric], reference_info, ...
            p.Results.GeometryTolerance);

        if is_noddi
            assert_same_geometry(metric_info, support_info, ...
                p.Results.GeometryTolerance, [subject ' ' metric ' versus support']);
            supported_values = double(metric_map(support_mask));
            if any(~isfinite(supported_values))
                error('generate_all_tract_metric_tables:NonFiniteNODDI', ...
                    '%s %s has non-finite values inside fitting support.', ...
                    subject, metric);
            end
            if any(supported_values < -p.Results.ValueTolerance | ...
                    supported_values > 1 + p.Results.ValueTolerance)
                error('generate_all_tract_metric_tables:NODDIOutOfBounds', ...
                    '%s %s has values outside [0,1] inside fitting support.', ...
                    subject, metric);
            end
            analysis_support = support_mask;
            zero_supported = nnz(supported_values == 0);
        else
            analysis_support = true(size(metric_map));
            zero_supported = NaN;
        end

        values = summarize_metric(metric_map, analysis_support, lesion_mask, ...
            masks, metric, startsWith(subject, 'sub-C'), mask_summary.Complete);
        metric_rows.(metric){end+1,1} = make_row(subject, columns, values);
        metric_support_voxels = support_voxels;
        if ~is_noddi
            metric_support_voxels = NaN;
        end
        qc_records(end+1) = make_qc_record(subject, metric, 'PASS', 1, ... %#ok<AGROW>
            mask_summary, map_affine_error, metric_support_voxels, lesion_status, ...
            zero_supported);
        clear metric_map supported_values analysis_support
    end
    clear masks lesion_mask support_mask support_data

    if mod(si, 10) == 0 || si == numel(subject_ids)
        fprintf('  processed %d/%d subjects\n', si, numel(subject_ids));
    end
end

tables = struct();
for mi = 1:numel(metrics)
    metric = metrics{mi};
    if isempty(metric_rows.(metric))
        tables.(metric) = empty_tract_table(columns);
    else
        tables.(metric) = vertcat(metric_rows.(metric){:});
    end
end

if has_noddi && ~explicit_subjects && p.Results.ExpectedNODDISubjects > 0
    for metric = {'NDI','ODI','FWF'}
        if ~ismember(metric{1}, metrics)
            continue;
        end
        observed = height(tables.(metric{1}));
        if observed ~= p.Results.ExpectedNODDISubjects
            error('generate_all_tract_metric_tables:UnexpectedNODDICohort', ...
                'Expected %d %s rows; observed %d.', ...
                p.Results.ExpectedNODDISubjects, metric{1}, observed);
        end
    end
end

validation_records = empty_validation_records();
for mi = 1:numel(metrics)
    metric = metrics{mi};
    if ismember(metric, {'T1','MTR','FA','MD'}) && ~isempty(legacy_stats_dir)
        reference_file = fullfile(legacy_stats_dir, ...
            ['GroupTract' metric '_All.xlsx']);
        if isfile(reference_file)
            validation_records(end+1) = compare_reference_table( ... %#ok<AGROW>
                tables.(metric), reference_file, metric, explicit_subjects, ...
                p.Results.ValidationAbsoluteTolerance, ...
                p.Results.ValidationRelativeTolerance);
        end
    elseif ismember(metric, {'NDI','ODI','FWF'}) && ...
            ~isempty(reference_noddi_dir)
        reference_file = fullfile(reference_noddi_dir, ...
            ['GroupTract' metric '_All.csv']);
        if ~isfile(reference_file)
            error('generate_all_tract_metric_tables:MissingReferenceTable', ...
                'Reference NODDI table not found: %s', reference_file);
        end
        validation_records(end+1) = compare_reference_table( ... %#ok<AGROW>
            tables.(metric), reference_file, metric, explicit_subjects, ...
            p.Results.ValidationAbsoluteTolerance, ...
            p.Results.ValidationRelativeTolerance);
    end
end
validation = struct2table(validation_records);
if ~isempty(validation) && any(~validation.Passed)
    failed = strjoin(cellstr(validation.Metric(~validation.Passed)), ', ');
    error('generate_all_tract_metric_tables:ReferenceValidationFailed', ...
        'Generated tables failed reference validation: %s.', failed);
end

files = table('Size', [0, 4], ...
    'VariableTypes', {'string','string','string','double'}, ...
    'VariableNames', {'Metric','Format','Path','Subjects'});
for mi = 1:numel(metrics)
    metric = metrics{mi};
    table_value = tables.(metric);
    csv_path = fullfile(output_dir, ['GroupTract' metric '_All.csv']);
    atomic_write_lossless_numeric_csv(table_value, csv_path, p.Results.Force);
    files = [files; {string(display_metric(metric)), "CSV", ...
        string(csv_path), height(table_value)}]; %#ok<AGROW>
    if ismember(metric, {'T1','MTR','FA','MD'})
        xlsx_path = fullfile(output_dir, ['GroupTract' metric '_All.xlsx']);
        atomic_write_verified_numeric_xlsx(table_value, xlsx_path, ...
            p.Results.Force);
        files = [files; {string(display_metric(metric)), "XLSX", ...
            string(xlsx_path), height(table_value)}]; %#ok<AGROW>
    end
end

qc = struct2table(qc_records);
if p.Results.WriteQC
    atomic_writetable(qc, fullfile(output_dir, ...
        'metric_table_generation_QC.csv'), p.Results.Force);
    atomic_writetable(validation, fullfile(output_dir, ...
        'metric_table_reference_validation.csv'), p.Results.Force);
end

fprintf('Generated %d aggregate table files in %s.\n', height(files), output_dir);
if ~isempty(validation)
    fprintf('Reference checks: %d/%d PASS.\n', ...
        nnz(validation.Passed), height(validation));
end
end


function metrics = canonical_metrics(values)
values = upper(cellstr(string(values)));
allowed = {'T1','MTR','FA','MD','NDI','ODI','ISOVF','FWF'};
if any(~ismember(values, allowed))
    invalid = values(~ismember(values, allowed));
    error('generate_all_tract_metric_tables:UnknownMetric', ...
        'Unknown metric(s): %s', strjoin(invalid, ', '));
end
values(strcmp(values, 'ISOVF')) = {'FWF'};
metrics = unique(values, 'stable');
end


function subjects = normalize_subjects(values)
subjects = cellstr(string(values));
subjects = cellfun(@strtrim, subjects, 'UniformOutput', false);
for ii = 1:numel(subjects)
    if ~startsWith(subjects{ii}, 'sub-')
        subjects{ii} = ['sub-' subjects{ii}];
    end
end
subjects = unique(subjects, 'stable');
subjects = sort(subjects);
end


function subjects = discover_subjects(source_root, legacy_stats_dir)
% The final manuscript workbook is the cohort authority when available. It
% prevents an otherwise valid but excluded BIDS directory from silently
% entering a rerun (for example, a later QC exclusion retained on disk).
cohort_file = fullfile(legacy_stats_dir, 'GroupTractFA_All.xlsx');
if isfile(cohort_file)
    cohort = readtable(cohort_file, 'VariableNamingRule', 'preserve');
    if ~ismember('SubjectID', cohort.Properties.VariableNames)
        error('generate_all_tract_metric_tables:MissingCohortSubjectID', ...
            '%s has no SubjectID column.', cohort_file);
    end
    subjects = normalize_subjects(string(cohort.SubjectID));
    if numel(subjects) ~= height(cohort)
        error('generate_all_tract_metric_tables:DuplicateCohortSubjectID', ...
            '%s contains duplicate SubjectIDs.', cohort_file);
    end
    return;
end
candidates = {
    fullfile(source_root, 'derivatives', 'TractoFlow_post'), ...
    fullfile(source_root, 'derivatives', 'maps')};
subjects = {};
for ci = 1:numel(candidates)
    if ~isfolder(candidates{ci})
        continue;
    end
    listing = dir(fullfile(candidates{ci}, 'sub-*'));
    listing = listing([listing.isdir]);
    if ~isempty(listing)
        subjects = sort({listing.name});
        return;
    end
end
end


function path_value = find_reference_path(source_root, subject)
candidates = {
    fullfile(source_root, 'derivatives', 'TractoFlow', 'ses-01', subject, ...
        'Resample_T1', [subject '__t1_resampled.nii.gz']), ...
    fullfile(source_root, 'derivatives', 'maps', subject, ...
        [subject '_space-individual_FA.nii.gz']), ...
    fullfile(source_root, 'derivatives', 'maps', subject, ...
        [subject '_ses-01_space-individual_T1map.nii.gz'])};
path_value = first_existing_file(candidates);
if isempty(path_value)
    error('generate_all_tract_metric_tables:MissingReference', ...
        'No T1-grid reference found for %s.', subject);
end
end


function path_value = metric_map_path(source_root, noddi_root, subject, metric)
switch metric
    case 'T1'
        path_value = fullfile(source_root, 'derivatives', 'maps', subject, ...
            [subject '_ses-01_space-individual_T1map.nii.gz']);
    case 'MTR'
        path_value = fullfile(source_root, 'derivatives', 'maps', subject, ...
            [subject '_ses-01_space-individual_MTRmap.nii.gz']);
    case {'FA','MD'}
        path_value = fullfile(source_root, 'derivatives', 'maps', subject, ...
            [subject '_space-individual_' metric '.nii.gz']);
    case {'NDI','ODI','FWF'}
        path_value = fullfile(noddi_root, subject, ...
            [subject '_space-individual_' metric '.nii.gz']);
    otherwise
        error('generate_all_tract_metric_tables:InternalMetric', ...
            'Unhandled metric: %s', metric);
end
end


function [masks, summary] = load_legacy_masks(source_root, subject, reference_info, tolerance)
groups = tract_groups();
hemispheres = {'L','R',''};
masks = struct();
missing = strings(0, 1);
max_affine_error = 0;
min_union_dice = 1;
left_x = nan(numel(groups), 1);
right_x = nan(numel(groups), 1);

for gi = 1:numel(groups)
    group = groups{gi};
    for hi = 1:numel(hemispheres)
        hemi = hemispheres{hi};
        key = mask_key(group, hemi);
        path_value = find_mask_path(source_root, subject, [group hemi]);
        if isempty(path_value)
            masks.(key) = [];
            missing(end+1,1) = string([group hemi]); %#ok<AGROW>
            continue;
        end
        info = read_header_checked(path_value, ...
            [subject ' ' group hemi ' tract mask'], tolerance);
        raw = niftiread(info);
        oriented = flip(flip(permute(raw, [1 3 2]), 2), 3) > 0;
        if ~isequal(size3(oriented), size3_from_info(reference_info))
            error('generate_all_tract_metric_tables:MaskShapeMismatch', ...
                '%s %s%s reorients to %s; reference is %s.', ...
                subject, group, hemi, mat2str(size3(oriented)), ...
                mat2str(size3_from_info(reference_info)));
        end
        % Keep an existing-but-empty mask: the original manuscript helper
        % emitted four zeros for legacy metrics. NODDI later rejects the
        % same empty supported ROI instead of manufacturing a measurement.
        effective_affine = legacy_oriented_affine(info);
        affine_error = max(abs(effective_affine - reference_info.Transform.T), [], 'all');
        if affine_error > tolerance
            error('generate_all_tract_metric_tables:MaskAffineMismatch', ...
                ['%s %s%s effective affine after the manuscript''s permute/flip ' ...
                 'differs from the T1 grid by %.9g (limit %.9g).'], ...
                subject, group, hemi, affine_error, tolerance);
        end
        max_affine_error = max(max_affine_error, affine_error);
        masks.(key) = oriented;
    end

    if ~isempty(masks.(mask_key(group, 'L'))) && ...
            ~isempty(masks.(mask_key(group, 'R'))) && ...
            ~isempty(masks.(mask_key(group, ''))) && ...
            any(masks.(mask_key(group, 'L')), 'all') && ...
            any(masks.(mask_key(group, 'R')), 'all') && ...
            any(masks.(mask_key(group, '')), 'all')
        left = masks.(mask_key(group, 'L'));
        right = masks.(mask_key(group, 'R'));
        combined = masks.(mask_key(group, ''));
        union_mask = left | right;
        denominator = nnz(union_mask) + nnz(combined);
        if denominator == 0
            union_dice = 1;
        else
            union_dice = 2 * nnz(union_mask & combined) / denominator;
        end
        if union_dice < 0.99
            error('generate_all_tract_metric_tables:MaskUnionMismatch', ...
                '%s %s L/R union versus combined-mask Dice is %.6f.', ...
                subject, group, union_dice);
        end
        min_union_dice = min(min_union_dice, union_dice);
        left_x(gi) = world_centroid_x(left, reference_info.Transform.T);
        right_x(gi) = world_centroid_x(right, reference_info.Transform.T);
        if ~(left_x(gi) < right_x(gi))
            error('generate_all_tract_metric_tables:LeftRightReversal', ...
                ['%s %s failed anatomical L/R centroid check after the legacy ' ...
                 'permute/flip: L=%g mm, R=%g mm.'], ...
                subject, group, left_x(gi), right_x(gi));
        end
        if strcmp(group, 'Association') && ...
                ~(left_x(gi) < 0 && right_x(gi) > 0)
            error('generate_all_tract_metric_tables:AssociationDoesNotStraddleMidline', ...
                ['%s Association masks do not straddle anatomical x=0 after ' ...
                 'the legacy permute/flip: L=%g mm, R=%g mm.'], ...
                subject, left_x(gi), right_x(gi));
        end
    end
end

summary = struct();
summary.Complete = isempty(missing);
summary.Missing = strjoin(cellstr(missing), ';');
summary.MaxAffineError = max_affine_error;
summary.MinUnionDice = min_union_dice;
summary.MaxLeftCentroidX = max(left_x, [], 'omitnan');
summary.MinRightCentroidX = min(right_x, [], 'omitnan');
end


function path_value = find_mask_path(source_root, subject, stem)
candidates = {
    fullfile(source_root, 'derivatives', 'TractoFlow_post', subject, ...
        'groupTract', [stem '.nii.gz']), ...
    fullfile(source_root, 'derivatives', 'TractoFlow_post', subject, ...
        [stem '.nii.gz'])};
path_value = first_existing_file(candidates);
end


function [lesion, status] = load_lesion_mask(source_root, subject, reference_info, strict, tolerance)
reference_size = size3_from_info(reference_info);
if startsWith(subject, 'sub-C')
    lesion = false(reference_size);
    status = 'NOT_APPLICABLE_CONTROL';
    return;
end
path_value = fullfile(source_root, 'derivatives', 'lesionMask', subject, ...
    'ses-01', [subject '_ses-01_desc-lesionManual_mask.nii.gz']);
if ~isfile(path_value)
    if strict
        error('generate_all_tract_metric_tables:MissingLesionMask', ...
            '%s manual lesion mask not found: %s', subject, path_value);
    end
    lesion = false(reference_size);
    status = 'MISSING_LEGACY_TREATED_AS_NO_LESION';
    return;
end
[data, ~] = read_volume_checked(path_value, [subject ' lesion mask'], ...
    reference_info, tolerance);
lesion = data > 0;
status = 'PASS';
end


function values = summarize_metric(metric_map, support, lesion, masks, metric, is_control, all_masks_complete)
groups = tract_groups();
hemispheres = {'L','R',''};
values = nan(1, numel(tract_columns()));
cursor = 1;
is_noddi = ismember(metric, {'NDI','ODI','FWF'});

for gi = 1:numel(groups)
    for hi = 1:numel(hemispheres)
        tract = masks.(mask_key(groups{gi}, hemispheres{hi}));
        if isempty(tract)
            if is_noddi || all_masks_complete
                error('generate_all_tract_metric_tables:UnexpectedMissingMask', ...
                    'A required tract mask is missing during %s extraction.', metric);
            end
            % Exact legacy behavior: a missing tract returns four zeros.
            values(cursor:cursor+3) = 0;
            cursor = cursor + 4;
            continue;
        end

        roi_all = tract & support;
        if is_noddi && ~any(roi_all(:))
            error('generate_all_tract_metric_tables:EmptySupportedTract', ...
                '%s has an empty supported tract ROI.', metric);
        end
        all_values = finite_values(metric_map, roi_all, is_noddi, metric);
        if isempty(all_values)
            if is_noddi
                error('generate_all_tract_metric_tables:EmptyNODDIValues', ...
                    '%s has no valid supported tract values.', metric);
            end
            values(cursor:cursor+3) = 0;
            cursor = cursor + 4;
            continue;
        end

        if is_control
            if is_noddi
                lesion_value = NaN;
            else
                lesion_value = 0;
            end
            nawm_value = mean(all_values);
        else
            lesion_values = finite_values(metric_map, roi_all & lesion, is_noddi, metric);
            nawm_values = finite_values(metric_map, roi_all & ~lesion, is_noddi, metric);
            if isempty(lesion_values)
                lesion_value = NaN;
                if ~is_noddi
                    lesion_value = 0;
                end
            else
                lesion_value = mean(lesion_values);
            end
            if isempty(nawm_values)
                if is_noddi
                    error('generate_all_tract_metric_tables:EmptyNODDINAWM', ...
                        '%s has an empty supported NAWM tract ROI.', metric);
                end
                nawm_value = 0;
            else
                nawm_value = mean(nawm_values);
            end
        end

        all_value = mean(all_values);
        tail_value = prctile(all_values, tail_percentile(metric));
        values(cursor:cursor+3) = [lesion_value, nawm_value, all_value, tail_value];
        cursor = cursor + 4;
    end
end
if is_noddi
    finite_summaries = values(isfinite(values));
    if any(finite_summaries < 0 | finite_summaries > 1)
        error('generate_all_tract_metric_tables:NODDIAggregateOutOfBounds', ...
            '%s produced an aggregate outside the strict [0,1] domain.', metric);
    end
end
end


function values = finite_values(metric_map, roi, require_all_finite, metric)
values = metric_map(roi);
if require_all_finite && any(~isfinite(values))
    error('generate_all_tract_metric_tables:NonFiniteROI', ...
        '%s contains non-finite values in a supported ROI.', metric);
end
values = values(isfinite(values));
if require_all_finite
    % NODDI zero is a valid biological/model value; conversion to double is
    % explicit and never coupled to missingness. Legacy metrics deliberately
    % retain their on-disk datatype because paper_gen_excel_stats.m did so.
    values = double(values);
end
end


function percentile = tail_percentile(metric)
switch metric
    case {'MTR','FA','NDI'}
        percentile = 10;
    case {'T1','MD','ODI','FWF'}
        percentile = 90;
    otherwise
        error('generate_all_tract_metric_tables:MissingTailRule', ...
            'No tail rule is defined for %s.', metric);
end
end


function T = make_row(subject, columns, values)
T = array2table(values, 'VariableNames', columns);
T = addvars(T, string(subject), 'Before', 1, 'NewVariableNames', 'SubjectID');
end


function T = empty_tract_table(columns)
T = array2table(nan(0, numel(columns)), 'VariableNames', columns);
T = addvars(T, strings(0,1), 'Before', 1, 'NewVariableNames', 'SubjectID');
end


function columns = tract_columns()
groups = tract_groups();
hemispheres = {'L','R',''};
summaries = {'Lesion','NAWM','All','Tail'};
columns = cell(1, numel(groups) * numel(hemispheres) * numel(summaries));
cursor = 1;
for gi = 1:numel(groups)
    for hi = 1:numel(hemispheres)
        for si = 1:numel(summaries)
            columns{cursor} = [groups{gi} hemispheres{hi} '_' summaries{si}];
            cursor = cursor + 1;
        end
    end
end
end


function groups = tract_groups()
groups = {'Association','Cerebellar','Occipitoparietal','ProjectionBrainstem'};
end


function key = mask_key(group, hemi)
if isempty(hemi)
    hemi = 'Combined';
end
key = [group '_' hemi];
end


function x = world_centroid_x(mask, transform)
indices = find(mask);
[i, j, k] = ind2sub(size(mask), indices);
centroid = [mean(i) mean(j) mean(k) 1] * transform;
x = centroid(1);
end


function effective = legacy_oriented_affine(info)
% Map zero-based oriented indices back to the raw mask indices before the
% original flip(flip(permute(mask,[1 3 2]),2),3) operation.
n = size3_from_info(info);
index_transform = [1 0 0 0; 0 0 -1 0; 0 -1 0 0; ...
    0 n(2)-1 n(3)-1 1];
effective = index_transform * info.Transform.T;
end


function [data, info, affine_error] = read_volume_checked(path_value, label, reference_info, tolerance)
if ~isfile(path_value)
    error('generate_all_tract_metric_tables:MissingFile', ...
        '%s not found: %s', label, path_value);
end
info = read_header_checked(path_value, label, tolerance);
affine_error = assert_same_geometry(info, reference_info, tolerance, label);
data = niftiread(info);
data_size = size(data);
if numel(data_size) > 3 && any(data_size(4:end) ~= 1)
    error('generate_all_tract_metric_tables:NotThreeDimensional', ...
        '%s is not a scalar 3D image.', label);
end
data = reshape(data, size3_from_info(info));
end


function info = read_header_checked(path_value, label, tolerance)
if ~isfile(path_value)
    error('generate_all_tract_metric_tables:MissingFile', ...
        '%s not found: %s', label, path_value);
end
info = niftiinfo(path_value);
if numel(info.ImageSize) < 3 || any(~isfinite(double(info.ImageSize(1:3))))
    error('generate_all_tract_metric_tables:InvalidHeader', ...
        '%s has invalid dimensions.', label);
end
transform = info.Transform.T;
if ~isequal(size(transform), [4 4]) || any(~isfinite(transform(:))) || ...
        abs(det(transform(1:3,1:3))) < eps
    error('generate_all_tract_metric_tables:InvalidAffine', ...
        '%s has a missing, non-finite, or singular affine.', label);
end
if isfield(info, 'raw')
    raw = info.raw;
    q_code = double(field_or(raw, 'qform_code', 0));
    s_code = double(field_or(raw, 'sform_code', 0));
    if q_code <= 0 && s_code <= 0
        error('generate_all_tract_metric_tables:MissingSpatialTransform', ...
            '%s has neither a valid qform nor sform.', label);
    end
    if q_code > 0 && s_code > 0
        qform = nifti_qform(raw);
        sform = nifti_sform(raw);
        q_s_error = max(abs(qform - sform), [], 'all');
        q_axes = orientation_signature(qform);
        s_axes = orientation_signature(sform);
        if ~strcmp(q_axes, s_axes) || ...
                sign(det(qform(1:3,1:3))) ~= sign(det(sform(1:3,1:3)))
            error('generate_all_tract_metric_tables:QformSformMismatch', ...
                ['%s qform and sform disagree in orientation/handedness ' ...
                 '(q=%s, s=%s).'], label, q_axes, s_axes);
        end
        % Some source lesion masks contain sub-voxel qform/sform numeric
        % differences while both encode the same LIA orientation.  MATLAB
        % and the original analysis select the sform in these files, so a
        % numeric-only difference is documented rather than rejected.
        if q_s_error > max(1e-3, 10 * tolerance)
            warning('generate_all_tract_metric_tables:QformSformNumericDifference', ...
                ['%s qform/sform differ numerically by %.6g mm but retain ' ...
                 'the same %s orientation; selected affine is unchanged.'], ...
                label, q_s_error, s_axes);
        end
    end
end
end


function error_value = assert_same_geometry(info, reference, tolerance, label)
if ~isequal(size3_from_info(info), size3_from_info(reference))
    error('generate_all_tract_metric_tables:ImageShapeMismatch', ...
        '%s shape %s differs from reference %s.', label, ...
        mat2str(size3_from_info(info)), mat2str(size3_from_info(reference)));
end
if any(abs(double(info.PixelDimensions(1:3)) - ...
        double(reference.PixelDimensions(1:3))) > tolerance)
    error('generate_all_tract_metric_tables:VoxelSizeMismatch', ...
        '%s voxel sizes differ from the reference.', label);
end
error_value = max(abs(info.Transform.T - reference.Transform.T), [], 'all');
if error_value > tolerance
    error('generate_all_tract_metric_tables:ImageAffineMismatch', ...
        '%s affine differs from reference by %.9g (limit %.9g).', ...
        label, error_value, tolerance);
end
if sign(det(info.Transform.T(1:3,1:3))) ~= ...
        sign(det(reference.Transform.T(1:3,1:3)))
    error('generate_all_tract_metric_tables:AffineHandednessMismatch', ...
        '%s affine handedness differs from the reference.', label);
end
end


function matrix = nifti_sform(raw)
matrix = [double(raw.srow_x(:)'); double(raw.srow_y(:)'); ...
    double(raw.srow_z(:)'); 0 0 0 1];
end


function matrix = nifti_qform(raw)
b = double(raw.quatern_b);
c = double(raw.quatern_c);
d = double(raw.quatern_d);
square_sum = b*b + c*c + d*d;
if square_sum > 1 + 1e-5
    error('generate_all_tract_metric_tables:InvalidQuaternion', ...
        'NIfTI qform quaternion has norm greater than one.');
end
a = sqrt(max(0, 1 - square_sum));
rotation = [ ...
    a*a+b*b-c*c-d*d, 2*b*c-2*a*d,     2*b*d+2*a*c; ...
    2*b*c+2*a*d,     a*a+c*c-b*b-d*d, 2*c*d-2*a*b; ...
    2*b*d-2*a*c,     2*c*d+2*a*b,     a*a+d*d-c*c-b*b];
pixdim = double(raw.pixdim);
qfac = pixdim(1);
if qfac == 0
    qfac = 1;
end
matrix = eye(4);
matrix(1:3,1:3) = rotation * diag([pixdim(2), pixdim(3), pixdim(4)*qfac]);
matrix(1:3,4) = [double(raw.qoffset_x); double(raw.qoffset_y); double(raw.qoffset_z)];
end


function signature = orientation_signature(matrix)
linear = matrix(1:3,1:3);
[~, world_axes] = max(abs(linear), [], 1);
if numel(unique(world_axes)) ~= 3
    error('generate_all_tract_metric_tables:ObliqueOrientationAmbiguous', ...
        'Could not derive a unique NIfTI orientation signature.');
end
letters = repmat('?', 1, 3);
for axis = 1:3
    direction = sign(linear(world_axes(axis), axis));
    if direction == 0
        error('generate_all_tract_metric_tables:ZeroOrientationAxis', ...
            'A NIfTI affine has a zero dominant orientation axis.');
    end
    switch world_axes(axis)
        case 1
            if direction < 0, letters(axis) = 'L'; else, letters(axis) = 'R'; end
        case 2
            if direction < 0, letters(axis) = 'P'; else, letters(axis) = 'A'; end
        case 3
            if direction < 0, letters(axis) = 'I'; else, letters(axis) = 'S'; end
    end
end
signature = letters;
end


function value = field_or(structure, field_name, default_value)
if isfield(structure, field_name)
    value = structure.(field_name);
else
    value = default_value;
end
end


function result = compare_reference_table(generated, reference_file, metric, allow_reference_superset, atol, rtol)
expected = readtable(reference_file, 'VariableNamingRule', 'preserve');
required = [{'SubjectID'}, tract_columns()];
if any(~ismember(required, expected.Properties.VariableNames))
    error('generate_all_tract_metric_tables:ReferenceSchemaMismatch', ...
        '%s does not contain the locked 49-column tract schema.', reference_file);
end
expected.SubjectID = string(expected.SubjectID);
generated.SubjectID = string(generated.SubjectID);
reference_subject_count = height(expected);
if numel(unique(expected.SubjectID)) ~= height(expected) || ...
        numel(unique(generated.SubjectID)) ~= height(generated)
    error('generate_all_tract_metric_tables:DuplicateSubject', ...
        '%s comparison contains duplicate SubjectIDs.', metric);
end

[present, positions] = ismember(generated.SubjectID, expected.SubjectID);
if any(~present)
    error('generate_all_tract_metric_tables:ReferenceSubjectMissing', ...
        '%s reference is missing %d generated subjects.', metric, nnz(~present));
end
if ~allow_reference_superset && ...
        ~isequal(sort(generated.SubjectID), sort(expected.SubjectID))
    error('generate_all_tract_metric_tables:ReferenceCohortMismatch', ...
        '%s generated and reference SubjectID sets differ.', metric);
end
expected = expected(positions, :);

% Historical workbooks may contain all-zero/all-missing placeholders for the
% two recovered grouped-mask profiles. Exclude those rows when validating an
% historical baseline; a final promoted baseline has no such migrations.
recovered_migration = false(height(generated), 1);
if ~allow_reference_superset && any(strcmp(metric, {'T1','FA','MD'}))
    names = tract_columns();
    target_values = double(expected{:, names});
    observed_values = double(generated{:, names});
    reference_placeholder = all(isnan(target_values) | target_values == 0, 2);
    observed_measurement = any(isfinite(observed_values) & observed_values ~= 0, 2);
    recovered_migration = reference_placeholder & observed_measurement;
    migration_count = nnz(recovered_migration);
    assert(migration_count == 0 || migration_count == 2, ...
        ['Expected either the final promoted baseline (zero migrations) or ' ...
         'the historical baseline (two regenerated grouped-mask rows).']);
end

compared = 0;
failed = 0;
max_error = 0;
missingness_failed = 0;
for column = tract_columns()
    observed = double(generated.(column{1}));
    target = double(expected.(column{1}));
    missing_mismatch = xor(isnan(observed), isnan(target)) & ...
        ~recovered_migration;
    missingness_failed = missingness_failed + nnz(missing_mismatch);
    finite = isfinite(observed) & isfinite(target);
    differences = abs(observed(finite) - target(finite));
    allowed = atol + rtol .* abs(target(finite));
    migrated_finite = recovered_migration(finite);
    audited_differences = differences(~migrated_finite);
    failed = failed + nnz((differences > allowed) & ~migrated_finite) + ...
        nnz(missing_mismatch) + ...
        nnz(xor(isinf(observed), isinf(target)) & ~recovered_migration);
    compared = compared + nnz(~recovered_migration);
    if ~isempty(audited_differences)
        max_error = max(max_error, max(audited_differences));
    end
end
result = struct('Metric', string(display_metric(metric)), ...
    'ReferenceFile', string(reference_file), ...
    'GeneratedSubjects', height(generated), ...
    'ReferenceSubjects', reference_subject_count, ...
    'ComparedCells', compared, ...
    'FailedCells', failed, ...
    'MissingnessFailures', missingness_failed, ...
    'MaximumAbsoluteError', max_error, ...
    'Passed', failed == 0);
end


function records = empty_qc_records()
prototype = struct('SubjectID', string.empty, 'Metric', string.empty, ...
    'Status', string.empty, 'RowWritten', NaN, 'MasksComplete', false, ...
    'MissingMasks', string.empty, 'MaximumMaskAffineError', NaN, ...
    'MinimumLRUnionDice', NaN, 'MaximumLeftCentroidX_mm', NaN, ...
    'MinimumRightCentroidX_mm', NaN, 'MapAffineError', NaN, ...
    'SupportVoxels', NaN, 'ZeroValuesInsideSupport', NaN, ...
    'LesionMaskStatus', string.empty);
records = repmat(prototype, 0, 1);
end


function record = make_qc_record(subject, metric, status, row_written, masks, map_error, support_voxels, lesion_status, zeros)
record = struct('SubjectID', string(subject), ...
    'Metric', string(display_metric(metric)), 'Status', string(status), ...
    'RowWritten', row_written, 'MasksComplete', masks.Complete, ...
    'MissingMasks', string(masks.Missing), ...
    'MaximumMaskAffineError', masks.MaxAffineError, ...
    'MinimumLRUnionDice', masks.MinUnionDice, ...
    'MaximumLeftCentroidX_mm', masks.MaxLeftCentroidX, ...
    'MinimumRightCentroidX_mm', masks.MinRightCentroidX, ...
    'MapAffineError', map_error, 'SupportVoxels', support_voxels, ...
    'ZeroValuesInsideSupport', zeros, ...
    'LesionMaskStatus', string(lesion_status));
end


function records = empty_validation_records()
prototype = struct('Metric', string.empty, 'ReferenceFile', string.empty, ...
    'GeneratedSubjects', NaN, 'ReferenceSubjects', NaN, ...
    'ComparedCells', NaN, 'FailedCells', NaN, ...
    'MissingnessFailures', NaN, 'MaximumAbsoluteError', NaN, ...
    'Passed', false);
records = repmat(prototype, 0, 1);
end


function atomic_writetable(T, final_path, force)
if isfile(final_path) && ~force
    error('generate_all_tract_metric_tables:OutputExists', ...
        'Output already exists (use Force=true to replace it): %s', final_path);
end
[output_dir, name, extension] = fileparts(final_path);
temporary = fullfile(output_dir, [name '.partial' extension]);
cleanup = onCleanup(@() delete_if_exists(temporary));
writetable(T, temporary);
if isfile(final_path)
    delete(final_path);
end
[moved, message] = movefile(temporary, final_path, 'f');
if ~moved
    error('generate_all_tract_metric_tables:AtomicMoveFailed', ...
        'Could not finalize %s: %s', final_path, message);
end
clear cleanup
end


function atomic_write_lossless_numeric_csv(T, final_path, force)
% Write model input without writetable's version-dependent float formatting.
validate_model_table_schema(T, final_path);
if isfile(final_path) && ~force
    error('generate_all_tract_metric_tables:OutputExists', ...
        'Output already exists (use Force=true to replace it): %s', final_path);
end
[output_dir, name, extension] = fileparts(final_path);
temporary = fullfile(output_dir, [name '.partial' extension]);
cleanup = onCleanup(@() delete_if_exists(temporary));

file_id = fopen(temporary, 'wt');
if file_id < 0
    error('generate_all_tract_metric_tables:CSVOpenFailed', ...
        'Could not open temporary CSV for writing: %s', temporary);
end
file_cleanup = onCleanup(@() close_if_open(file_id));
fprintf(file_id, '%s\n', strjoin(T.Properties.VariableNames, ','));
numeric_values = T{:, 2:end};
subject_ids = string(T{:, 1});
for row = 1:height(T)
    subject = char(subject_ids(row));
    if isempty(subject) || ~isempty(regexp(subject, '[,\"\r\n]', 'once'))
        error('generate_all_tract_metric_tables:UnsafeSubjectID', ...
            'SubjectID cannot be written unambiguously to CSV: %s', subject);
    end
    fprintf(file_id, '%s', subject);
    for column = 1:size(numeric_values, 2)
        value = numeric_values(row, column);
        if isinf(value)
            error('generate_all_tract_metric_tables:InfiniteModelValue', ...
                '%s contains an infinite value at row %d, column %d.', ...
                final_path, row, column + 1);
        elseif isnan(value)
            fprintf(file_id, ',NaN');
        else
            fprintf(file_id, ',%.17g', value);
        end
    end
    fprintf(file_id, '\n');
end
clear file_cleanup

assert_numeric_table_readback_identity(T, temporary, final_path);
finalize_temporary_file(temporary, final_path);
clear cleanup
end


function atomic_write_verified_numeric_xlsx(T, final_path, force)
% XLSX remains the established legacy input; explicitly audit its round trip.
validate_model_table_schema(T, final_path);
if isfile(final_path) && ~force
    error('generate_all_tract_metric_tables:OutputExists', ...
        'Output already exists (use Force=true to replace it): %s', final_path);
end
[output_dir, name, extension] = fileparts(final_path);
temporary = fullfile(output_dir, [name '.partial' extension]);
cleanup = onCleanup(@() delete_if_exists(temporary));
writetable(T, temporary);

observed = readtable(temporary, 'VariableNamingRule', 'preserve', ...
    'TextType', 'string');
assert_same_table_schema_and_ids(T, observed, final_path);
expected_values = double(T{:, 2:end});
observed_values = double(observed{:, 2:end});
if ~isequaln(expected_values, observed_values)
    error('generate_all_tract_metric_tables:XLSXRoundTripMismatch', ...
        'Exact XLSX round-trip identity failed for %s.', final_path);
end

finalize_temporary_file(temporary, final_path);
clear cleanup
end


function validate_model_table_schema(T, path_value)
if width(T) < 2 || ~strcmp(T.Properties.VariableNames{1}, 'SubjectID')
    error('generate_all_tract_metric_tables:ModelTableSchema', ...
        '%s must have SubjectID first followed by numeric columns.', path_value);
end
if ~(isstring(T{:,1}) || iscellstr(T{:,1}) || ischar(T{:,1}))
    error('generate_all_tract_metric_tables:ModelTableSubjectType', ...
        '%s SubjectID must be text.', path_value);
end
for column = 2:width(T)
    if ~isnumeric(T{:, column})
        error('generate_all_tract_metric_tables:ModelTableNumericType', ...
            '%s column %s is not numeric.', path_value, ...
            T.Properties.VariableNames{column});
    end
end
end


function assert_numeric_table_readback_identity(expected, csv_path, label)
observed = readtable(csv_path, 'VariableNamingRule', 'preserve', ...
    'TextType', 'string');
assert_same_table_schema_and_ids(expected, observed, label);
if ~isequaln(double(expected{:, 2:end}), double(observed{:, 2:end}))
    error('generate_all_tract_metric_tables:CSVRoundTripMismatch', ...
        'Lossless CSV read-back identity failed for %s.', label);
end
end


function assert_same_table_schema_and_ids(expected, observed, label)
if ~isequal(expected.Properties.VariableNames, observed.Properties.VariableNames) || ...
        height(expected) ~= height(observed)
    error('generate_all_tract_metric_tables:SerializedSchemaMismatch', ...
        'Serialized schema or row count differs for %s.', label);
end
if ~isequal(string(expected{:,1}), string(observed{:,1}))
    error('generate_all_tract_metric_tables:SerializedSubjectMismatch', ...
        'Serialized SubjectID order differs for %s.', label);
end
end


function finalize_temporary_file(temporary, final_path)
if isfile(final_path)
    delete(final_path);
end
[moved, message] = movefile(temporary, final_path, 'f');
if ~moved
    error('generate_all_tract_metric_tables:AtomicMoveFailed', ...
        'Could not finalize %s: %s', final_path, message);
end
end


function close_if_open(file_id)
if file_id >= 0
    fclose(file_id);
end
end


function delete_if_exists(path_value)
if isfile(path_value)
    delete(path_value);
end
end


function value = display_metric(metric)
if strcmp(metric, 'FWF')
    value = 'ISOVF';
else
    value = metric;
end
end


function value = first_existing_file(candidates)
value = '';
for ii = 1:numel(candidates)
    if isfile(candidates{ii})
        value = candidates{ii};
        return;
    end
end
end


function value = require_dir(path_value, label)
if isempty(path_value) || ~isfolder(path_value)
    error('generate_all_tract_metric_tables:MissingDirectory', ...
        '%s directory not found: %s', label, path_value);
end
value = path_value;
end


function result = is_text_scalar(value)
result = ischar(value) || (isstring(value) && isscalar(value));
end


function result = size3(value)
result = size(value);
result(end+1:3) = 1;
result = result(1:3);
end


function result = size3_from_info(info)
result = double(info.ImageSize);
result(end+1:3) = 1;
result = result(1:3);
end
