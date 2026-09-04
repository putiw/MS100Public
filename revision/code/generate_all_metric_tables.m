function outputs = generate_all_metric_tables(varargin)
% GENERATE_ALL_METRIC_TABLES One MATLAB entry point for all paper qMRI tables.
%
% The AMICO fit and ANTs resampling must already be complete. This function
% performs only the manuscript's MATLAB aggregation step and writes no images
% or participant-level maps.
%
% Minimal full-cohort call:
%   generate_all_metric_tables( ...
%       'SourceRoot', '/read/only/MsBIDS_or_snapshot', ...
%       'NODDIMapsDir', '/local/work/transformed_maps', ...
%       'OutputDir', '/local/work/metric_tables');
%
% That call generates:
%   * GroupTract T1/MTR/FA/MD (XLSX + CSV)
%   * GroupTract NDI/ODI/FWF (CSV; FWF is reported as ISOVF)
%   * ClassicalRegionNODDI_All.csv
%   * ClassicalRegionAllMetrics.csv (all seven metrics)
%
% Legacy tract and classical values are actively recomputed and must match
% the original manuscript workbooks. ReferenceNODDIDir is optional and is
% intended only for the one-time Python-to-MATLAB migration comparison.

p = inputParser;
addParameter(p, 'SourceRoot', '', @is_text_scalar);
addParameter(p, 'NODDIMapsDir', '', @is_text_scalar);
addParameter(p, 'OutputDir', '', @is_text_scalar);
addParameter(p, 'LegacyStatsDir', '', @is_text_scalar);
addParameter(p, 'ReferenceNODDIDir', '', @is_text_scalar);
addParameter(p, 'FreesurferDir', '', @is_text_scalar);
addParameter(p, 'FreesurferMatlabDir', '', @is_text_scalar);
addParameter(p, 'Metrics', {'T1','MTR','FA','MD','NDI','ODI','ISOVF'}, ...
    @(x) iscellstr(x) || isstring(x));
addParameter(p, 'Subjects', {}, @(x) iscellstr(x) || isstring(x));
addParameter(p, 'IncludeClassical', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'GeometryTolerance', 1e-4, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'ExpectedSourceSubjects', 132, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ExpectedNODDISubjects', 132, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'Force', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'WriteQC', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

source_root = require_dir(char(p.Results.SourceRoot), 'SourceRoot');
noddi_maps_dir = require_dir(char(p.Results.NODDIMapsDir), 'NODDIMapsDir');
output_dir = char(p.Results.OutputDir);
if isempty(output_dir)
    error('generate_all_metric_tables:MissingOutputDir', 'OutputDir is required.');
end

script_dir = fileparts(mfilename('fullpath'));
package_root = fileparts(script_dir);
helper_dir = fullfile(package_root, 'helper', 'matlab');
addpath(helper_dir, '-begin');
assert_safe_output(output_dir, source_root, package_root);

legacy_stats_dir = char(p.Results.LegacyStatsDir);
if isempty(legacy_stats_dir)
    legacy_stats_dir = fullfile(source_root, 'derivatives', 'derivativesStats');
end
require_dir(legacy_stats_dir, 'LegacyStatsDir');

[tract_files, tract_qc, tract_validation] = ...
    generate_all_tract_metric_tables( ...
        'SourceRoot', source_root, ...
        'NODDIMapsDir', noddi_maps_dir, ...
        'OutputDir', output_dir, ...
        'Metrics', p.Results.Metrics, ...
        'Subjects', p.Results.Subjects, ...
        'LegacyStatsDir', legacy_stats_dir, ...
        'ReferenceNODDIDir', p.Results.ReferenceNODDIDir, ...
        'GeometryTolerance', p.Results.GeometryTolerance, ...
        'ExpectedSourceSubjects', p.Results.ExpectedSourceSubjects, ...
        'ExpectedNODDISubjects', p.Results.ExpectedNODDISubjects, ...
        'Force', p.Results.Force, ...
        'WriteQC', p.Results.WriteQC);

outputs = struct();
outputs.TractFiles = tract_files;
outputs.TractQC = tract_qc;
outputs.TractValidation = tract_validation;

if p.Results.IncludeClassical
    freesurfer_dir = char(p.Results.FreesurferDir);
    if isempty(freesurfer_dir)
        freesurfer_dir = fullfile(source_root, 'derivatives', 'freesurfer');
    end
    require_dir(freesurfer_dir, 'FreesurferDir');
    legacy_classical = fullfile(legacy_stats_dir, 'icometrixLesionStats.xlsx');
    require_file(legacy_classical, 'legacy classical-region table');
    noddi_output = fullfile(output_dir, 'ClassicalRegionNODDI_All.csv');
    combined_output = fullfile(output_dir, 'ClassicalRegionAllMetrics.csv');

    [classical_noddi, classical_qc, classical_validation, classical_all] = ...
        gen_excel_stats_ico_noddi( ...
            'FreesurferMatlabDir', p.Results.FreesurferMatlabDir, ...
            'FreesurferDir', freesurfer_dir, ...
            'SourceSnapshot', source_root, ...
            'NODDIMapsDir', noddi_maps_dir, ...
            'OutputFile', noddi_output, ...
            'CombinedOutputFile', combined_output, ...
            'LegacyTable', legacy_classical, ...
            'Subjects', p.Results.Subjects, ...
            'RequireFullCohort', isempty(p.Results.Subjects), ...
            'GeometryTolerance', p.Results.GeometryTolerance, ...
            'Force', p.Results.Force, ...
            'WriteQC', p.Results.WriteQC);
    outputs.ClassicalNODDI = classical_noddi;
    outputs.ClassicalAllMetrics = classical_all;
    outputs.ClassicalQC = classical_qc;
    outputs.ClassicalValidation = classical_validation;
else
    outputs.ClassicalNODDI = table();
    outputs.ClassicalAllMetrics = table();
    outputs.ClassicalQC = table();
    outputs.ClassicalValidation = table();
end

fprintf('\nAll requested metric tables passed validation.\n');
fprintf('Output directory: %s\n', output_dir);
end


function assert_safe_output(output_dir, source_root, package_root)
output = canonical_path(output_dir);
source = canonical_path(source_root);
package = canonical_path(package_root);
if strcmp(output, '/') || strcmp(output, char(java.lang.System.getProperty('user.home')))
    error('generate_all_metric_tables:UnsafeOutput', ...
        'OutputDir cannot be a filesystem or home-directory root.');
end
if strcmp(output, '/Volumes') || startsWith(output, ['/Volumes' filesep])
    error('generate_all_metric_tables:OutputOnSourceVolume', ...
        'OutputDir must be local and cannot be below /Volumes.');
end
if is_same_or_below(output, source) || is_same_or_below(source, output)
    error('generate_all_metric_tables:SourceOutputOverlap', ...
        'SourceRoot and OutputDir must be separate trees.');
end
if is_same_or_below(output, package) || is_same_or_below(package, output)
    error('generate_all_metric_tables:PackageOutputOverlap', ...
        'OutputDir must be outside the public reproduction package.');
end
end


function result = is_same_or_below(path_value, parent)
result = strcmp(path_value, parent) || ...
    startsWith(path_value, [parent filesep]);
end


function value = canonical_path(path_value)
value = char(java.io.File(path_value).getCanonicalPath());
end


function path_value = require_dir(path_value, label)
if isempty(path_value) || ~isfolder(path_value)
    error('generate_all_metric_tables:MissingDirectory', ...
        '%s directory not found: %s', label, path_value);
end
end


function require_file(path_value, label)
if ~isfile(path_value)
    error('generate_all_metric_tables:MissingFile', ...
        '%s not found: %s', label, path_value);
end
end


function result = is_text_scalar(value)
result = ischar(value) || (isstring(value) && isscalar(value));
end
