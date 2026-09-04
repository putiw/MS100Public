function varargout = noddi_reproduction_support(action, varargin)
%NODDI_REPRODUCTION_SUPPORT Internal validation/provenance support.
% This private helper keeps reviewer safeguards separate from the
% rev_qMRI.m-derived analysis so that the scientific code remains easy to
% compare with its public source. It is not a user entry point.

switch action
    case 'add_msfc_sdmt'
        [varargout{1}, varargout{2}] = add_msfc_sdmt(varargin{1});
    case 'input_provenance'
        varargout{1} = input_provenance(varargin{1}, varargin{2});
    case 'classical_input_provenance'
        varargout{1} = classical_input_provenance(varargin{1}, varargin{2});
    case 'validate_input_cohorts'
        varargout{1} = validate_input_cohorts(varargin{1}, varargin{2});
    case 'compare_legacy_rows'
        varargout{1} = compare_legacy_rows(varargin{1});
    case 'validate_recovered_mask_migration'
        varargout{1} = validate_recovered_mask_migration(varargin{1});
    case 'write_json'
        write_json(varargin{1}, varargin{2});
    otherwise
        error('Unknown NODDI reproduction-support action: %s', action);
end
end

function audit = classical_input_provenance(classicalPath, legacyPath)
% Bind the new aggregate to the exact legacy classical cohort without
% exposing participant identifiers in a public report.
require_file(classicalPath);
require_file(legacyPath);
classical = readtable(classicalPath);
legacy = readtable(legacyPath);
assert(ismember('SubjectID',classical.Properties.VariableNames) && ...
       ismember('SubjectID',legacy.Properties.VariableNames), ...
    'Classical/legacy input is missing SubjectID.');
classicalIds = string(classical.SubjectID);
legacyIds = string(legacy.SubjectID);
exactIds = isequal(sort(classicalIds),sort(legacyIds));
assert(height(classical) == 87 && numel(unique(classicalIds)) == 87, ...
    'Classical NODDI aggregate is not the 87-subject unique cohort.');
assert(height(legacy) == 87 && numel(unique(legacyIds)) == 87, ...
    'Legacy classical table is not the 87-subject unique cohort.');
assert(exactIds, ...
    'Classical NODDI SubjectIDs do not exactly match the legacy classical cohort.');
details = dir(classicalPath);
audit = struct('name',get_filename(classicalPath), ...
    'sha256',sha256_file(classicalPath), 'size_bytes',details.bytes, ...
    'rows',height(classical), 'unique_subjects',numel(unique(classicalIds)), ...
    'legacy_rows',height(legacy), 'exact_legacy_subject_set',exactIds, ...
    'subject_ids_disclosed',false, 'status','PASS');
end

function result = compare_legacy_rows(actualRows)
actual = struct2table(actualRows);
reference = historical_reference_rows();
records = struct([]);
for ri = 1:height(reference)
    sheet = char(reference.Sheet(ri));
    adjusted = strcmp(sheet, 'Demo');
    outcome = char(reference.Outcome(ri));
    metric = char(reference.Metric(ri));
    ai = strcmp(actual.Outcome, outcome) & strcmp(actual.Metric, metric) & ...
        actual.Adjusted == adjusted;
    assert(nnz(ai) == 1, 'Could not identify one %s/%s/%s row.', ...
        sheet, outcome, metric);
    if ismember(outcome, {'EDSS','MSPro'})
        observedPerformance = sprintf('%.3f [%.3f, %.3f]', ...
            actual.Performance(ai), actual.CI_lower(ai), actual.CI_upper(ai));
    else
        observedPerformance = sprintf('%.3f', actual.Performance(ai));
    end
    observedAic = sprintf('%.2f', actual.AIC(ai));
    record = struct('Sheet',sheet, 'Outcome',outcome, 'Metric',metric, ...
        'N_expected',reference.N(ri), 'N_observed',actual.N(ai), ...
        'p_expected',reference.p(ri), 'p_observed',actual.p_raw(ai), ...
        'p_abs_difference',abs(reference.p(ri)-actual.p_raw(ai)), ...
        'Performance_expected',char(reference.Performance(ri)), ...
        'Performance_observed',observedPerformance, ...
        'AIC_expected',char(reference.AIC(ri)), 'AIC_observed',observedAic);
    record.p_display_expected = format_p_for_audit(reference.p(ri));
    record.p_display_observed = format_p_for_audit(actual.p_raw(ai));
    record.All_match = record.N_expected == record.N_observed && ...
        strcmp(record.p_display_expected, record.p_display_observed) && ...
        (reference.p(ri) < 0.05) == (actual.p_raw(ai) < 0.05) && ...
        strcmp(strtrim(record.Performance_expected), ...
            strtrim(record.Performance_observed)) && ...
        strcmp(strtrim(record.AIC_expected), strtrim(record.AIC_observed));
    records = append_struct(records, record);
end
result = struct();
result.historical_validation_results = [ ...
    'embedded 20-row FA/MD manuscript-display reference in ' ...
    'helper/matlab/private/noddi_reproduction_support.m'];
result.historical_validation_source_workbook_sha256 = ...
    '9d6b66a2288c2f6bd5c6e046c9b3213ef18d09dffc678689d253ad9f3fe21432';
result.analysis_source_code = 'MS100Public/rev_qMRI.m';
result.analysis_source_code_sha256 = ...
    'df4feae5033fbed00b5d749140560be9ddc516333f9552da1700d4abf3bd8c45';
result.rows_checked = numel(records);
result.rows_matched = sum([records.All_match]);
result.all_match = all([records.All_match]);
result.maximum_absolute_p_difference = max([records.p_abs_difference]);
result.match_definition = [ ...
    'N and significance exact; displayed AUC/R2, AIC, and p exact; ' ...
    'full-precision p difference recorded diagnostically'];
result.records = records;
end


function value = format_p_for_audit(pValue)
if pValue < 0.0001
    value = '<0.0001';
else
    value = sprintf('%.4f', pValue);
end
end


function reference = historical_reference_rows()
sheet = [repmat("Results",10,1); repmat("Demo",10,1)];
outcome = repmat(["EDSS";"EDSS";"MSPro";"MSPro";"T25FW";"T25FW"; ...
    "9HPT-D";"9HPT-D";"9HPT-ND";"9HPT-ND"], 2, 1);
metric = repmat(["FA";"MD"], 10, 1);
n = [79;79;78;78;125;125;128;128;128;128; ...
     79;79;78;78;125;125;128;128;128;128];
p = [ ...
    5.0284578255999e-6; 2.8822607631336652e-6; ...
    1.1767152332662559e-5; 3.6577542105429211e-6; ...
    8.2167683768119559e-10; 2.8310687127941489e-14; ...
    3.6400660263780082e-11; 2.7736516416965169e-9; ...
    3.33066907387547e-16; 9.1551096962168064e-23; ...
    4.707831045248908e-5; 9.1379579691409818e-5; ...
    9.0660068473165324e-4; 8.2712045329415632e-4; ...
    2.564101123647688e-7; 3.6093350530563839e-11; ...
    2.377725349944626e-9; 3.2663818216871482e-7; ...
    2.775557561562891e-14; 6.3282712403633923e-15];
performance = [ ...
    "0.897 [0.823, 0.959]"; "0.902 [0.826, 0.963]"; ...
    "0.861 [0.772, 0.931]"; "0.886 [0.803, 0.952]"; ...
    "0.265"; "0.376"; "0.295"; "0.245"; "0.412"; "0.434"; ...
    "0.935 [0.876, 0.978]"; "0.923 [0.859, 0.973]"; ...
    "0.966 [0.926, 0.993]"; "0.965 [0.917, 0.996]"; ...
    "0.343"; "0.431"; "0.322"; "0.267"; "0.456"; "0.469"];
aic = ["67.22";"66.40";"76.13";"68.11";"77.49";"56.94"; ...
    "-57.90";"-49.27";"-32.51";"-37.40"; ...
    "60.81";"67.11";"45.71";"47.06";"67.48";"49.58"; ...
    "-58.96";"-49.07";"-38.62";"-41.62"];
reference = table(sheet,outcome,metric,n,p,performance,aic, ...
    'VariableNames',{'Sheet','Outcome','Metric','N','p','Performance','AIC'});
end


function provenance = input_provenance(statsDir, metricsDir)
helperRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
configuration = jsondecode(fileread(fullfile( ...
    helperRoot, 'config', 'analysis_config.json')));
specifications = configuration.frozen_aggregate_inputs;
template = struct('name','', 'sha256','', 'expected_sha256','', ...
    'observed_fingerprint','', 'match_mode','', 'semantic_sha256','', ...
    'matches_completed_run',false, 'size_bytes',0);
provenance = repmat(template, numel(specifications), 1);
for i = 1:numel(specifications)
    if iscell(specifications)
        specification = specifications{i};
    else
        specification = specifications(i);
    end
    if strcmp(specification.directory, 'stats')
        path = fullfile(statsDir, specification.name);
        generatedPath = fullfile(metricsDir, specification.name);
        if startsWith(specification.name, 'GroupTract') && isfile(generatedPath)
            % Mirror the unified model: regenerated MATLAB tract workbooks
            % supersede their frozen manuscript copies when both are present.
            path = generatedPath;
        end
    else
        path = fullfile(metricsDir, specification.name);
    end
    require_file(path);
    details = dir(path);
    rawDigest = sha256_file(path);
    semanticDigest = '';
    if strcmp(specification.match_mode, 'raw_sha256')
        observed = rawDigest;
    elseif strcmp(specification.match_mode, 'semantic_legacy_tract_5sig_v2')
        semanticDigest = semantic_legacy_tract_5sig_sha256(path);
        observed = semanticDigest;
    elseif strcmp(specification.match_mode, 'semantic_noddi_tract_5dp_v1')
        semanticDigest = semantic_noddi_tract_5dp_sha256(path);
        observed = semanticDigest;
    else
        error('Unsupported input match mode: %s', specification.match_mode);
    end
    provenance(i) = struct('name',specification.name, 'sha256',rawDigest, ...
        'expected_sha256',specification.expected_fingerprint, ...
        'observed_fingerprint',observed, ...
        'match_mode',specification.match_mode, ...
        'semantic_sha256',semanticDigest, ...
        'matches_completed_run',strcmp(observed,specification.expected_fingerprint), ...
        'size_bytes',details.bytes);
end
end


function digest = semantic_noddi_tract_5dp_sha256(path)
% Canonicalize float32-derived NODDI aggregates at five decimal places.
% MATLAB and the historical Python extractor can differ by about 1e-8 in
% percentile interpolation even when they read the same maps. Five-decimal
% canonicalization is finer than manuscript display precision while
% still binding every value, missing cell, subject, and column to one digest.
cells = readcell(path);
assert(size(cells,2) == 49, ...
    'NODDI tract CSV must contain SubjectID plus 48 metric columns.');
nColumns = size(cells,2);
headers = strings(1,nColumns);
for column = 1:nColumns
    value = cells{1,column};
    assert(ischar(value) || (isstring(value) && isscalar(value)), ...
        'NODDI tract CSV has a non-text header.');
    headers(column) = strtrim(string(value));
end
assert(headers(1) == "SubjectID" && all(strlength(headers) > 0), ...
    'Invalid NODDI tract CSV header.');
assert(numel(unique(headers)) == nColumns, ...
    'Duplicate NODDI tract CSV headers.');

canonicalRows = strings(0,1);
subjects = strings(0,1);
for row = 2:size(cells,1)
    values = cells(row,1:nColumns);
    if all(cellfun(@semantic_missing, values))
        continue;
    end
    subjectValue = values{1};
    assert(ischar(subjectValue) || ...
        (isstring(subjectValue) && isscalar(subjectValue)), ...
        'NODDI tract row has no text SubjectID.');
    subject = strtrim(string(subjectValue));
    assert(~ismissing(subject) && strlength(subject) > 0, ...
        'NODDI tract row has no SubjectID.');
    assert(~any(subjects == subject), ...
        'Duplicate SubjectID in NODDI tract CSV.');
    subjects(end+1,1) = subject; %#ok<AGROW>
    encoded = strings(1,nColumns-1);
    for column = 2:nColumns
        value = values{column};
        if semantic_missing(value)
            encoded(column-1) = "NA";
        else
            assert(isnumeric(value) && isscalar(value) && isfinite(value), ...
                'Non-numeric or non-finite populated NODDI tract value.');
            numericValue = double(value);
            % Map-level QC alone permits a <=1e-6 numerical excursion.
            % Aggregates are never clipped and must satisfy the exact domain.
            assert(numericValue >= 0 && numericValue <= 1, ...
                'NODDI tract aggregate lies outside the strict [0,1] domain.');
            encoded(column-1) = string(sprintf('%.5f', numericValue));
        end
    end
    canonicalRows(end+1,1) = subject + sprintf('\t') + ...
        strjoin(encoded, sprintf('\t')); %#ok<AGROW>
end
helperRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
configuration = jsondecode(fileread(fullfile( ...
    helperRoot, 'config', 'analysis_config.json')));
expectedRows = double(configuration.study.expected_primary_tract_subjects);
assert(numel(canonicalRows) == expectedRows, ...
    'NODDI tract CSV must contain %d unique subject rows.', expectedRows);
canonicalRows = sort(canonicalRows);
payload = ['semantic_noddi_tract_5dp_v1' newline ...
    char(strjoin(headers, sprintf('\t'))) newline ...
    char(strjoin(canonicalRows, newline)) newline];
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(unicode2native(payload, 'UTF-8'));
raw = typecast(engine.digest(), 'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end


function result = validate_recovered_mask_migration(validation)
% The historical FA/MD model gate predates recovery of two grouped masks.
% Binary rows must remain exact because they use joint classical cases;
% continuous tract rows must gain exactly the two recovered participants.
records = validation.records;
outcomes = string({records.Outcome});
binary = ismember(outcomes, ["EDSS","MSPro"]);
continuous = ismember(outcomes, ["T25FW","9HPT-D","9HPT-ND"]);
assert(nnz(binary) == 8 && nnz(continuous) == 12 && ...
    all(binary | continuous), ...
    'Unexpected historical FA/MD validation row set.');
assert(all([records(binary).All_match]), ...
    'A binary FA/MD result changed during recovered-mask migration.');
expectedN = [records(continuous).N_expected];
observedN = [records(continuous).N_observed];
assert(all(observedN == expectedN + 2), ...
    'Continuous FA/MD rows did not gain exactly two recovered participants.');
result = struct('status','PASS', ...
    'interpretation','expected recovered-mask cohort migration', ...
    'binary_rows_unchanged',nnz(binary), ...
    'continuous_rows_with_two_added_subjects',nnz(continuous), ...
    'unexpected_rows',0);
end


function digest = semantic_legacy_tract_5sig_sha256(path)
% Canonicalize legacy aggregates at five significant figures. The maps are
% single precision, so a fresh MATLAB aggregation can legitimately differ
% from an old workbook by one float32 step. The table generator separately
% enforces the tighter cellwise 1e-10 + 1e-6 relative reference tolerance.
cells = readcell(path, 'Sheet', 1);
while ~isempty(cells) && all(cellfun(@semantic_missing, cells(:,end)))
    cells(:,end) = [];
end
assert(size(cells,2) >= 2, 'Legacy tract workbook has no metric columns.');
nColumns = size(cells,2);
headers = strings(1,nColumns);
for column = 1:nColumns
    value = cells{1,column};
    assert(ischar(value) || (isstring(value) && isscalar(value)), ...
        'Legacy tract workbook has a non-text header.');
    headers(column) = strtrim(string(value));
end
assert(headers(1) == "SubjectID" && all(strlength(headers) > 0), ...
    'Invalid legacy tract workbook header.');
assert(numel(unique(headers)) == nColumns, ...
    'Duplicate legacy tract workbook headers.');

canonicalRows = strings(0,1);
subjects = strings(0,1);
for row = 2:size(cells,1)
    values = cells(row,1:nColumns);
    if all(cellfun(@semantic_missing, values))
        continue;
    end
    measurements = values(2:end);
    if all(cellfun(@semantic_missing, measurements))
        % Normalize an identifier-only row whose entire metric profile is
        % blank; this is layout metadata, not an analyzed observation.
        continue;
    end
    subjectValue = values{1};
    assert(ischar(subjectValue) || (isstring(subjectValue) && isscalar(subjectValue)), ...
        'Populated legacy tract row has no text SubjectID.');
    subject = strtrim(string(subjectValue));
    assert(~ismissing(subject) && strlength(subject) > 0, ...
        'Populated legacy tract row has no SubjectID.');
    assert(~any(subjects == subject), 'Duplicate SubjectID in legacy tract workbook.');
    subjects(end+1,1) = subject; %#ok<AGROW>
    encoded = strings(1,nColumns-1);
    for column = 2:nColumns
        value = values{column};
        if semantic_missing(value)
            encoded(column-1) = "NA";
        else
            assert(isnumeric(value) && isscalar(value) && isfinite(value), ...
                'Non-numeric or non-finite populated legacy tract value.');
            encoded(column-1) = string(sprintf('%.5g', double(value)));
        end
    end
    canonicalRows(end+1,1) = subject + sprintf('\t') + ...
        strjoin(encoded, sprintf('\t')); %#ok<AGROW>
end
canonicalRows = sort(canonicalRows);
payload = ['semantic_legacy_tract_5sig_v2' newline ...
    char(strjoin(headers, sprintf('\t'))) newline ...
    char(strjoin(canonicalRows, newline)) newline];
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(unicode2native(payload, 'UTF-8'));
raw = typecast(engine.digest(), 'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end


function value = semantic_missing(item)
if isempty(item)
    value = true;
elseif ischar(item)
    value = isempty(strtrim(item));
elseif isstring(item)
    value = ismissing(item) || strlength(strtrim(item)) == 0;
elseif isnumeric(item)
    value = isscalar(item) && isnan(item);
else
    try
        value = ismissing(item);
        value = isscalar(value) && value;
    catch
        value = false;
    end
end
end


function audit = validate_input_cohorts(clinical, metricsDir)
helperRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
configuration = jsondecode(fileread(fullfile( ...
    helperRoot, 'config', 'analysis_config.json')));
study = configuration.study;
clinicalSubjects = string(clinical.SubjectID);
expectedClinical = double(study.expected_clinical_subjects);
expectedNoddi = double(study.expected_primary_tract_subjects);
expectedOverlap = double(study.expected_clinical_noddi_overlap);
expectedNoddiOnly = double(study.expected_noddi_only_subjects);
expectedClinicalOnly = double(study.expected_clinical_only_subjects);
assert(height(clinical) == expectedClinical && ...
    numel(unique(clinicalSubjects)) == expectedClinical, ...
    'Clinical input does not match the frozen unique-subject count.');
metricNames = {'NDI','ODI','FWF'};
referenceSubjects = strings(0,1);
for index = 1:numel(metricNames)
    metric = metricNames{index};
    tablePath = fullfile(metricsDir, ['GroupTract' metric '_All.csv']);
    require_file(tablePath);
    metricTable = readtable(tablePath);
    metricSubjects = string(metricTable.SubjectID);
    assert(height(metricTable) == expectedNoddi && ...
        numel(unique(metricSubjects)) == expectedNoddi, ...
        '%s input does not match the frozen unique-subject count.', metric);
    if index == 1
        referenceSubjects = sort(metricSubjects);
    else
        assert(isequal(sort(metricSubjects), referenceSubjects), ...
            'NDI, ODI, and FWF subject sets must match exactly.');
    end
end
overlap = numel(intersect(clinicalSubjects, referenceSubjects));
noddiOnly = numel(setdiff(referenceSubjects, clinicalSubjects));
clinicalOnly = numel(setdiff(clinicalSubjects, referenceSubjects));
assert(overlap == expectedOverlap && noddiOnly == expectedNoddiOnly && ...
    clinicalOnly == expectedClinicalOnly, ...
    'Clinical/NODDI subject overlap does not match the frozen cohort accounting.');
audit = struct('clinical_unique_subjects',expectedClinical, ...
    'noddi_unique_subjects',expectedNoddi, ...
    'clinical_noddi_overlap',overlap, ...
    'noddi_only_subjects',noddiOnly, ...
    'clinical_only_subjects',clinicalOnly, ...
    'excluded_for_absent_grouped_masks', ...
        double(study.expected_missing_original_group_mask_subjects), ...
    'metric_subject_sets_match',true, 'status','PASS');
end


function digest = sha256_file(path)
fid = fopen(path, 'rb');
assert(fid >= 0, 'Could not hash required input.');
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(bytes);
raw = typecast(engine.digest(), 'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));
end


function [T, audit] = add_msfc_sdmt(T)
required = {'SubjectID','Group','T25FW','x9HPTD','x9HPTND','SDMTcorrect'};
missing = setdiff(required, T.Properties.VariableNames);
assert(isempty(missing), 'MSFC source columns missing: %s', ...
    strjoin(missing, ', '));
subjects = string(T.SubjectID);
assert(numel(unique(subjects)) == height(T), ...
    'MSFC scoring requires one row per unique subject.');

components = [T.T25FW,T.x9HPTD,T.x9HPTND,T.SDMTcorrect];
valid = all(isfinite(components), 2) & ...
    all(components(:,1:3) > 0, 2) & components(:,4) >= 0;
reference = valid & ~startsWith(subjects, 'sub-C');
assert(height(T) == 132, 'Expected 132 clinical rows; found %d.', height(T));
assert(sum(valid) == 127, 'Expected 127 complete scores; found %d.', sum(valid));
assert(sum(reference) == 84, 'Expected 84 complete MS reference rows.');

arm = nan(height(T),1);
arm(valid) = (1./T.x9HPTD(valid) + 1./T.x9HPTND(valid))/2;
armMean = mean(arm(reference));
armSd = std(arm(reference));
legMean = mean(T.T25FW(reference));
legSd = std(T.T25FW(reference));
cogMean = mean(T.SDMTcorrect(reference));
cogSd = std(T.SDMTcorrect(reference));

zArm = nan(height(T),1);
zLeg = nan(height(T),1);
zCog = nan(height(T),1);
zArm(valid) = (arm(valid)-armMean)/armSd;
zLeg(valid) = -(T.T25FW(valid)-legMean)/legSd;
zCog(valid) = (T.SDMTcorrect(valid)-cogMean)/cogSd;
score = mean([zArm,zLeg,zCog],2,'omitmissing');
score(~valid) = NaN;
T.MSFC_SDMT = score;

audit = struct();
audit.endpoint = 'MSFC_SDMT';
audit.label = 'MSFC-SDMT (not the original PASAT-based MSFC)';
audit.direction = 'higher is better';
audit.reference_population = 'complete-case MS participants';
audit.reference_n = sum(reference);
audit.complete_score_n = sum(valid);
audit.total_n = height(T);
audit.missing_score_n = sum(~valid);
audit.motor_input_policy = ...
    'clinician-selected recorded fields used as stored; no reconstruction';
audit.missing_policy = 'complete case; no imputation';
audit.reference = struct( ...
    'arm_reciprocal_mean',armMean, 'arm_reciprocal_sd',armSd, ...
    'T25FW_seconds_mean',legMean, 'T25FW_seconds_sd',legSd, ...
    'SDMT_correct_mean',cogMean, 'SDMT_correct_sd',cogSd);
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


function write_json(path, value)
fid = fopen(path, 'w');
assert(fid >= 0, 'Could not open %s.', path);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(value, PrettyPrint=true));
end
