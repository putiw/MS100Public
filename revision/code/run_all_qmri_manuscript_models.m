function report = run_all_qmri_manuscript_models(varargin)
%RUN_ALL_QMRI_MANUSCRIPT_MODELS Refit every manuscript qMRI model from tables.
%
% This is the human-facing MATLAB entry point for the final combined
% analysis. It actively refits:
%   seven imaging metrics: T1, MTR, FA, MD, NDI, ODI, and FWF/ISOVF;
%   seven outcomes: EDSS, MSPro, T25FW, 9HPT-D, 9HPT-ND, SDMT, MSFC-SDMT;
%   tract-based and classical-region models; and
%   unadjusted and demographically adjusted specifications.
%
% Required name-value arguments:
%   statsDir      clinicalScore.xlsx and icometrixLesionStats.xlsx; original
%                 GroupTract workbooks are accepted as a compatibility input
%   metricsDir    actively regenerated GroupTract tables for all seven
%                 metrics plus ClassicalRegionAllMetrics.csv
%   classicalFile ClassicalRegionNODDI_All.csv compatibility view
%   outputDir     destination for model and figure-input CSVs
%
% Example:
%   run_all_qmri_manuscript_models( ...
%       'statsDir', '/path/to/derivativesStats', ...
%       'metricsDir', '/path/to/noddi/tables', ...
%       'classicalFile', '/path/to/ClassicalRegionNODDI_All.csv', ...
%       'outputDir', '/path/to/output');

script_dir = fileparts(mfilename('fullpath'));
package_root = fileparts(script_dir);
helper_dir = fullfile(package_root, 'helper', 'matlab');
addpath(helper_dir, '-begin');

report = rev_qMRI_noddi(varargin{:}, 'analysisScope', 'all-metrics');
end
