function [results, qc, validation, combined_results] = gen_excel_stats_ico_noddi(varargin)
% GEN_EXCEL_STATS_ICO_NODDI Generate classical-region qMRI lesion means.
%
% This is a minimal extension of the study's original classical-region
% MATLAB logic. That source is not a standalone file in the current public
% repository; the definitions below are explicit and validated against the
% original 87-row aggregate table.
% It actively recomputes the original T1/MTR/FA/MD values as a hard identity
% gate, then adds NDI/ODI/FWF (manuscript label: ISOVF) using the same four
% FreeSurfer-derived regions. NODDI values are restricted to valid AMICO
% fitting support, so a biological zero is retained while an empty region is
% stored as NaN. Voxel-map QC permits at most a 1e-6 numerical excursion from
% [0,1] and never clips; every reported aggregate must lie strictly in [0,1].
%
% Required name-value inputs:
%   FreesurferDir  - root containing sub-*/mri/{aparc+aseg,ribbon}.mgz
%   SourceSnapshot - BIDS root/snapshot containing lesion masks and maps
%   NODDIMapsDir   - root containing sub-* transformed NODDI maps
%   OutputFile     - destination ClassicalRegionNODDI_All.csv
%   LegacyTable    - icometrixLesionStats.xlsx; fixes cohort and row order
%
% Optional:
%   CombinedOutputFile - all seven metrics in one CSV; defaults to
%                        ClassicalRegionAllMetrics.csv beside OutputFile
%   FreesurferMatlabDir - FreeSurfer MATLAB utilities containing MRIread.m
%   Subjects          - explicit subset used for reproducibility tests
%   RequireFullCohort - require the manuscript's 87 subjects (default true)
%   GeometryTolerance - vox2ras0 tolerance in mm (default 1e-4)
%   Tolerance          - legacy value-validation tolerance (default 1e-12)
%   Force              - replace existing outputs (default false)
%   WriteQC            - write compact QC/validation CSVs (default true)
%
% The two model-facing CSVs use 17 significant digits and must survive an
% exact read-back identity check before either output is finalized.

p = inputParser;
addParameter(p, 'FreesurferMatlabDir', '', @is_text_scalar);
addParameter(p, 'FreesurferDir', '', @is_text_scalar);
addParameter(p, 'SourceSnapshot', '', @is_text_scalar);
addParameter(p, 'NODDIMapsDir', '', @is_text_scalar);
addParameter(p, 'OutputFile', '', @is_text_scalar);
addParameter(p, 'CombinedOutputFile', '', @is_text_scalar);
addParameter(p, 'LegacyTable', '', @is_text_scalar);
addParameter(p, 'Subjects', {}, @(x) iscellstr(x) || isstring(x));
addParameter(p, 'RequireFullCohort', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'GeometryTolerance', 1e-4, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'Tolerance', 1e-12, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ValueTolerance', 1e-6, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 0 && x <= 1e-6);
addParameter(p, 'Force', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'WriteQC', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

freesurfer_matlab_dir = char(p.Results.FreesurferMatlabDir);
standard_freesurfer_matlab_dir = '/Applications/freesurfer/7.4.1/matlab';
if isempty(freesurfer_matlab_dir) && isfolder(standard_freesurfer_matlab_dir)
    freesurfer_matlab_dir = standard_freesurfer_matlab_dir;
end
if ~isempty(freesurfer_matlab_dir)
    require_dir(freesurfer_matlab_dir, 'FreesurferMatlabDir');
    addpath(freesurfer_matlab_dir, '-begin');
end
if exist('MRIread', 'file') ~= 2
    error('gen_excel_stats_ico_noddi:MissingMRIread', ...
        ['FreeSurfer MRIread.m was not found. Supply FreesurferMatlabDir ' ...
         'pointing to the FreeSurfer MATLAB utilities directory.']);
end

freesurfer_dir = require_dir(char(p.Results.FreesurferDir), 'FreesurferDir');
snapshot_dir = require_dir(char(p.Results.SourceSnapshot), 'SourceSnapshot');
noddi_maps_dir = require_dir(char(p.Results.NODDIMapsDir), 'NODDIMapsDir');
output_file = char(p.Results.OutputFile);
legacy_file = char(p.Results.LegacyTable);
if isempty(output_file)
    error('gen_excel_stats_ico_noddi:MissingOutput', 'OutputFile is required.');
end
require_file(legacy_file, 'LegacyTable');
[output_dir, output_name] = fileparts(output_file);
combined_output_file = char(p.Results.CombinedOutputFile);
if isempty(combined_output_file)
    combined_output_file = fullfile(output_dir, 'ClassicalRegionAllMetrics.csv');
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

legacy_validation_file = fullfile(output_dir, ...
    [output_name '_Legacy_T1_MTR_FA_MD_validation.csv']);
compatibility_validation_file = fullfile(output_dir, ...
    [output_name '_Legacy_FA_MD_validation.csv']);
qc_file = fullfile(output_dir, [output_name '_QC.csv']);
planned = {output_file, combined_output_file};
if p.Results.WriteQC
    planned = [planned, {qc_file, legacy_validation_file, ...
        compatibility_validation_file}];
end
if ~p.Results.Force
    existing = planned(cellfun(@isfile, planned));
    if ~isempty(existing)
        error('gen_excel_stats_ico_noddi:OutputExists', ...
            'Refusing to replace existing output (use Force=true): %s', existing{1});
    end
end

regions = {'periventricular', 'juxtacortical', ...
    'infratentorial', 'deepwhitematter'};
legacy_metrics = {'T1','MTR','FA','MD'};
noddi_metrics = {'NDI','ODI','FWF'};
all_metrics = [legacy_metrics, noddi_metrics];
results = initialize_metric_table(regions, noddi_metrics);
combined_results = initialize_metric_table(regions, all_metrics);
qc = table('Size', [0, 9], ...
    'VariableTypes', {'string','string','double','double','double','logical', ...
        'double','double','string'}, ...
    'VariableNames', {'SubjectID','Region','LesionVoxels', ...
        'SupportedLesionVoxels','SupportFraction','HasSupportedLesion', ...
        'MaximumLegacyAffineError','MaximumNODDIAffineError','Status'});
validation = table('Size', [0, 8], ...
    'VariableTypes', {'string','string','string','double','double','double','logical','double'}, ...
    'VariableNames', {'SubjectID','Metric','Region','Calculated','Expected', ...
        'AbsoluteError','Passed','Tolerance'});

legacy = readtable(legacy_file, 'VariableNamingRule', 'preserve');
if ~ismember('SubjectID', legacy.Properties.VariableNames)
    error('gen_excel_stats_ico_noddi:MissingSubjectID', ...
        'LegacyTable must contain a SubjectID column.');
end
legacy.SubjectID = string(legacy.SubjectID);
if numel(unique(legacy.SubjectID)) ~= height(legacy)
    error('gen_excel_stats_ico_noddi:DuplicateSubjectID', ...
        'LegacyTable contains duplicate SubjectIDs.');
end
if p.Results.RequireFullCohort && height(legacy) ~= 87
    error('gen_excel_stats_ico_noddi:UnexpectedLegacyCohort', ...
        'Expected 87 manuscript subjects; observed %d.', height(legacy));
end
if any(startsWith(legacy.SubjectID, 'sub-C'))
    error('gen_excel_stats_ico_noddi:UnexpectedControl', ...
        'The legacy classical cohort unexpectedly contains a control SubjectID.');
end

if isempty(p.Results.Subjects)
    subject_ids = legacy.SubjectID;
else
    subject_ids = string(p.Results.Subjects(:));
    subject_ids(~startsWith(subject_ids, 'sub-')) = ...
        "sub-" + subject_ids(~startsWith(subject_ids, 'sub-'));
    subject_ids = unique(subject_ids, 'stable');
    if p.Results.RequireFullCohort && numel(subject_ids) ~= 87
        error('gen_excel_stats_ico_noddi:SubsetWithFullCohortGate', ...
            'Set RequireFullCohort=false when using an explicit subset.');
    end
    if any(~ismember(subject_ids, legacy.SubjectID))
        error('gen_excel_stats_ico_noddi:UnknownSubsetSubject', ...
            'Subjects contains an ID absent from LegacyTable.');
    end
end

for si = 1:numel(subject_ids)
    subject = char(subject_ids(si));
    legacy_index = find(legacy.SubjectID == string(subject), 1);
    fprintf('Processing classical qMRI subject: %s\n', subject);

    lesion_path = fullfile(snapshot_dir, 'derivatives', 'lesionMask', ...
        subject, 'ses-01', [subject '_ses-01_desc-lesionManual_mask.nii.gz']);
    aparc_path = fullfile(freesurfer_dir, subject, 'mri', 'aparc+aseg.mgz');
    ribbon_path = fullfile(freesurfer_dir, subject, 'mri', 'ribbon.mgz');
    support_path = fullfile(noddi_maps_dir, subject, ...
        [subject '_space-individual_NODDI-support.nii.gz']);
    require_file(lesion_path, [subject ' lesion mask']);
    require_file(aparc_path, [subject ' aparc+aseg']);
    require_file(ribbon_path, [subject ' ribbon']);
    require_file(support_path, [subject ' NODDI support']);

    lesion_vol = MRIread(lesion_path);
    aparc_vol = MRIread(aparc_path);
    ribbon_vol = MRIread(ribbon_path);
    support_vol = MRIread(support_path);
    lesion_mask = lesion_vol.vol > 0;
    support_mask = support_vol.vol > 0;
    assert_same_geometry(subject, lesion_vol, aparc_vol, ...
        p.Results.GeometryTolerance, 'aparc+aseg');
    assert_same_geometry(subject, lesion_vol, ribbon_vol, ...
        p.Results.GeometryTolerance, 'ribbon');
    assert_same_geometry(subject, lesion_vol, support_vol, ...
        p.Results.GeometryTolerance, 'NODDI support');
    if ~any(support_mask(:))
        error('gen_excel_stats_ico_noddi:EmptySupport', ...
            '%s has an empty transformed AMICO fitting support.', subject);
    end

    % These definitions preserve the study's original classical-region logic.
    region_masks = struct();
    region_masks.periventricular = ismember(aparc_vol.vol, [4, 43]);
    region_masks.infratentorial = ismember(aparc_vol.vol, [16, 7, 8]);
    region_masks.juxtacortical = ismember(ribbon_vol.vol, [3, 42]);
    region_masks.deepwhitematter = lesion_mask & ...
        ~region_masks.periventricular & ...
        ~region_masks.infratentorial & ...
        ~region_masks.juxtacortical;

    metric_maps = struct();
    maximum_legacy_affine_error = 0;
    maximum_noddi_affine_error = 0;
    for mi = 1:numel(all_metrics)
        metric = all_metrics{mi};
        map_path = classical_metric_path(snapshot_dir, noddi_maps_dir, subject, metric);
        require_file(map_path, [subject ' ' metric ' map']);
        metric_maps.(metric) = MRIread(map_path);
        affine_error = assert_same_geometry(subject, lesion_vol, ...
            metric_maps.(metric), p.Results.GeometryTolerance, [metric ' map']);
        if ismember(metric, legacy_metrics)
            maximum_legacy_affine_error = max(maximum_legacy_affine_error, affine_error);
        else
            maximum_noddi_affine_error = max(maximum_noddi_affine_error, affine_error);
            supported_values = double(metric_maps.(metric).vol(support_mask));
            if any(~isfinite(supported_values)) || ...
                    any(supported_values < -p.Results.ValueTolerance | ...
                        supported_values > 1 + p.Results.ValueTolerance)
                error('gen_excel_stats_ico_noddi:InvalidNODDIValues', ...
                    '%s %s contains non-finite or out-of-range supported values.', ...
                    subject, metric);
            end
        end
    end

    noddi_row = table(string(subject), 'VariableNames', {'SubjectID'});
    combined_row = table(string(subject), 'VariableNames', {'SubjectID'});
    for ri = 1:numel(regions)
        region = regions{ri};
        original_roi = region_masks.(region) & lesion_mask;
        noddi_roi = original_roi & support_mask;
        original_n = nnz(original_roi);
        supported_n = nnz(noddi_roi);
        fraction = supported_n / max(original_n, 1);
        has_values = supported_n > 0;
        status = "PASS";
        if ~has_values
            status = "EMPTY_SUPPORTED_LESION";
        end
        qc = [qc; {string(subject), string(region), original_n, supported_n, ...
            fraction, has_values, maximum_legacy_affine_error, ...
            maximum_noddi_affine_error, status}]; %#ok<AGROW>

        for mi = 1:numel(legacy_metrics)
            metric = legacy_metrics{mi};
            calculated = mean_or_nan(double(metric_maps.(metric).vol(original_roi)));
            combined_row.([region metric]) = calculated;
            expected = legacy{legacy_index, [region metric]};
            [passed, absolute_error] = compare_value(calculated, expected, ...
                p.Results.Tolerance);
            validation = [validation; {string(subject), string(metric), ...
                string(region), calculated, expected, absolute_error, ...
                passed, p.Results.Tolerance}]; %#ok<AGROW>
        end

        for mi = 1:numel(noddi_metrics)
            metric = noddi_metrics{mi};
            value = mean_or_nan(double(metric_maps.(metric).vol(noddi_roi)));
            if isfinite(value) && (value < 0 || value > 1)
                error('gen_excel_stats_ico_noddi:NODDIAggregateOutOfBounds', ...
                    ['%s %s %s produced an aggregate outside the strict ' ...
                     '[0,1] domain.'], subject, metric, region);
            end
            noddi_row.([region metric]) = value;
            combined_row.([region metric]) = value;
        end
    end
    results = [results; noddi_row]; %#ok<AGROW>
    combined_results = [combined_results; combined_row]; %#ok<AGROW>
end

if ~isempty(validation) && ~all(validation.Passed)
    failed = nnz(~validation.Passed);
    error('gen_excel_stats_ico_noddi:LegacyValidationFailed', ...
        'Legacy T1/MTR/FA/MD validation failed in %d/%d cells.', ...
        failed, height(validation));
end

atomic_write_lossless_numeric_csv(results, output_file, p.Results.Force);
atomic_write_lossless_numeric_csv(combined_results, combined_output_file, ...
    p.Results.Force);
if p.Results.WriteQC
    atomic_writetable(qc, qc_file, p.Results.Force);
    atomic_writetable(validation, legacy_validation_file, p.Results.Force);
    compatibility_validation = validation(ismember(validation.Metric, ["FA","MD"]), :);
    atomic_writetable(compatibility_validation, compatibility_validation_file, ...
        p.Results.Force);
end

fprintf('Classical all-metric rows: %d\n', height(combined_results));
fprintf('Supported NODDI region rows: %d/%d\n', ...
    nnz(qc.HasSupportedLesion), height(qc));
fprintf('Legacy T1/MTR/FA/MD validation: %d/%d PASS\n', ...
    nnz(validation.Passed), height(validation));
fprintf('NODDI table: %s\n', output_file);
fprintf('Combined seven-metric table: %s\n', combined_output_file);
end


function T = initialize_metric_table(regions, metrics)
num_cols = 1 + numel(regions) * numel(metrics);
types = [{'string'}, repmat({'double'}, 1, num_cols - 1)];
names = {'SubjectID'};
for ri = 1:numel(regions)
    for mi = 1:numel(metrics)
        names{end+1} = [regions{ri} metrics{mi}]; %#ok<AGROW>
    end
end
T = table('Size', [0, num_cols], 'VariableTypes', types, ...
    'VariableNames', names);
end


function path_value = classical_metric_path(snapshot, noddi_root, subject, metric)
switch metric
    case 'T1'
        path_value = fullfile(snapshot, 'derivatives', 'maps', subject, ...
            [subject '_ses-01_space-individual_T1map.nii.gz']);
    case 'MTR'
        path_value = fullfile(snapshot, 'derivatives', 'maps', subject, ...
            [subject '_ses-01_space-individual_MTRmap.nii.gz']);
    case {'FA','MD'}
        path_value = fullfile(snapshot, 'derivatives', 'maps', subject, ...
            [subject '_space-individual_' metric '.nii.gz']);
    case {'NDI','ODI','FWF'}
        path_value = fullfile(noddi_root, subject, ...
            [subject '_space-individual_' metric '.nii.gz']);
    otherwise
        error('gen_excel_stats_ico_noddi:UnknownMetric', ...
            'Unknown classical-region metric: %s', metric);
end
end


function affine_error = assert_same_geometry(subject, reference, observed, tolerance, label)
if ~isequal(size(reference.vol), size(observed.vol))
    error('gen_excel_stats_ico_noddi:GeometryMismatch', ...
        '%s %s size %s differs from reference %s.', subject, label, ...
        mat2str(size(observed.vol)), mat2str(size(reference.vol)));
end
if ~isfield(reference, 'vox2ras0') || ~isfield(observed, 'vox2ras0') || ...
        ~isequal(size(reference.vox2ras0), [4 4]) || ...
        ~isequal(size(observed.vox2ras0), [4 4])
    error('gen_excel_stats_ico_noddi:MissingAffine', ...
        '%s %s lacks a usable vox2ras0 affine.', subject, label);
end
affine_error = max(abs(double(reference.vox2ras0) - ...
    double(observed.vox2ras0)), [], 'all');
if affine_error > tolerance
    error('gen_excel_stats_ico_noddi:AffineMismatch', ...
        '%s %s affine difference is %.9g (limit %.9g).', ...
        subject, label, affine_error, tolerance);
end
if sign(det(reference.vox2ras0(1:3,1:3))) ~= ...
        sign(det(observed.vox2ras0(1:3,1:3)))
    error('gen_excel_stats_ico_noddi:AffineHandednessMismatch', ...
        '%s %s affine handedness differs from the lesion reference.', ...
        subject, label);
end
end


function value = mean_or_nan(values)
if isempty(values)
    value = NaN;
    return;
end
if any(~isfinite(values))
    values = values(isfinite(values));
end
if isempty(values)
    value = NaN;
else
    value = mean(values);
end
end


function [passed, absolute_error] = compare_value(calculated, expected, tolerance)
if isnan(calculated) && isnan(expected)
    absolute_error = 0;
    passed = true;
elseif xor(isnan(calculated), isnan(expected))
    absolute_error = Inf;
    passed = false;
else
    absolute_error = abs(calculated - expected);
    passed = absolute_error <= tolerance;
end
end


function atomic_writetable(T, final_path, force)
if isfile(final_path) && ~force
    error('gen_excel_stats_ico_noddi:OutputExists', ...
        'Output already exists (use Force=true): %s', final_path);
end
[output_dir, name, extension] = fileparts(final_path);
if ~isfolder(output_dir)
    mkdir(output_dir);
end
temporary = fullfile(output_dir, [name '.partial' extension]);
cleanup = onCleanup(@() delete_if_exists(temporary));
writetable(T, temporary);
if isfile(final_path)
    delete(final_path);
end
[moved, message] = movefile(temporary, final_path, 'f');
if ~moved
    error('gen_excel_stats_ico_noddi:AtomicMoveFailed', ...
        'Could not finalize %s: %s', final_path, message);
end
clear cleanup
end


function atomic_write_lossless_numeric_csv(T, final_path, force)
validate_model_table_schema(T, final_path);
if isfile(final_path) && ~force
    error('gen_excel_stats_ico_noddi:OutputExists', ...
        'Output already exists (use Force=true): %s', final_path);
end
[output_dir, name, extension] = fileparts(final_path);
if ~isfolder(output_dir)
    mkdir(output_dir);
end
temporary = fullfile(output_dir, [name '.partial' extension]);
cleanup = onCleanup(@() delete_if_exists(temporary));

file_id = fopen(temporary, 'wt');
if file_id < 0
    error('gen_excel_stats_ico_noddi:CSVOpenFailed', ...
        'Could not open temporary CSV for writing: %s', temporary);
end
file_cleanup = onCleanup(@() close_if_open(file_id));
fprintf(file_id, '%s\n', strjoin(T.Properties.VariableNames, ','));
numeric_values = T{:, 2:end};
subject_ids = string(T{:, 1});
for row = 1:height(T)
    subject = char(subject_ids(row));
    if isempty(subject) || ~isempty(regexp(subject, '[,\"\r\n]', 'once'))
        error('gen_excel_stats_ico_noddi:UnsafeSubjectID', ...
            'SubjectID cannot be written unambiguously to CSV: %s', subject);
    end
    fprintf(file_id, '%s', subject);
    for column = 1:size(numeric_values, 2)
        value = numeric_values(row, column);
        if isinf(value)
            error('gen_excel_stats_ico_noddi:InfiniteModelValue', ...
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

observed = readtable(temporary, 'VariableNamingRule', 'preserve', ...
    'TextType', 'string');
if ~isequal(T.Properties.VariableNames, observed.Properties.VariableNames) || ...
        height(T) ~= height(observed)
    error('gen_excel_stats_ico_noddi:SerializedSchemaMismatch', ...
        'Serialized schema or row count differs for %s.', final_path);
end
if ~isequal(string(T{:,1}), string(observed{:,1}))
    error('gen_excel_stats_ico_noddi:SerializedSubjectMismatch', ...
        'Serialized SubjectID order differs for %s.', final_path);
end
if ~isequaln(double(T{:, 2:end}), double(observed{:, 2:end}))
    error('gen_excel_stats_ico_noddi:CSVRoundTripMismatch', ...
        'Lossless CSV read-back identity failed for %s.', final_path);
end

if isfile(final_path)
    delete(final_path);
end
[moved, message] = movefile(temporary, final_path, 'f');
if ~moved
    error('gen_excel_stats_ico_noddi:AtomicMoveFailed', ...
        'Could not finalize %s: %s', final_path, message);
end
clear cleanup
end


function validate_model_table_schema(T, path_value)
if width(T) < 2 || ~strcmp(T.Properties.VariableNames{1}, 'SubjectID')
    error('gen_excel_stats_ico_noddi:ModelTableSchema', ...
        '%s must have SubjectID first followed by numeric columns.', path_value);
end
if ~(isstring(T{:,1}) || iscellstr(T{:,1}) || ischar(T{:,1}))
    error('gen_excel_stats_ico_noddi:ModelTableSubjectType', ...
        '%s SubjectID must be text.', path_value);
end
for column = 2:width(T)
    if ~isnumeric(T{:, column})
        error('gen_excel_stats_ico_noddi:ModelTableNumericType', ...
            '%s column %s is not numeric.', path_value, ...
            T.Properties.VariableNames{column});
    end
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


function path_value = require_dir(path_value, label)
if isempty(path_value) || ~isfolder(path_value)
    error('gen_excel_stats_ico_noddi:MissingDirectory', ...
        '%s directory not found: %s', label, path_value);
end
end


function require_file(path_value, label)
if ~isfile(path_value)
    error('gen_excel_stats_ico_noddi:MissingFile', ...
        '%s not found: %s', label, path_value);
end
end


function result = is_text_scalar(value)
result = ischar(value) || (isstring(value) && isscalar(value));
end
