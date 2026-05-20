function paper_gen_excel_stats(varargin)
% PAPER_GEN_EXCEL_STATS Generates statistics for various metrics in lesions, NAWM, and all tract
%
% This version includes measurements for:
%   - Lesion: voxels within lesions in the tract
%   - NAWM: voxels within tract but outside lesions
%   - All: all voxels within the tract (lesion + NAWM combined)
%
% Usage:
%   paper_gen_excel_stats('T1')              % Generate T1 statistics using manual lesion mask and dwi tracts
%   paper_gen_excel_stats('FA')              % Generate FA statistics using manual lesion mask and dwi tracts
%   paper_gen_excel_stats('MD')              % Generate MD statistics using manual lesion mask and dwi tracts
%   paper_gen_excel_stats('MTR')             % Generate MTR statistics using manual lesion mask and dwi tracts
%   paper_gen_excel_stats('AD')              % Generate AD statistics using manual lesion mask and dwi tracts
%   paper_gen_excel_stats('RD')              % Generate RD statistics using manual lesion mask and dwi tracts
%   paper_gen_excel_stats('T1', 'template')  % Use template tract masks
%   paper_gen_excel_stats('T1', 'lesionMask', 'LSTAI') % Use LSTAI lesion mask
%   paper_gen_excel_stats('T1', 'template', 'lesionMask', 'LSTAI') % Use both template and LSTAI
%
% Input Data:
%   - Lesion masks: derivatives/lesionMask/sub-*/ses-01/sub-*_ses-01_desc-lesionManual_mask.nii.gz
%   - Lesion masks (LSTAI): derivatives/lesionMask/sub-*/ses-01/sub-*_ses-01_desc-LSTAI_mask.nii.gz
%   - Metric maps: derivatives/maps/sub-*/sub-*_space-individual_*.nii.gz

    % Parse input arguments
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));
    p = inputParser;
    addRequired(p, 'metric', @(x) ismember(x, {'T1', 'FA', 'MD', 'MTR', 'AD', 'RD'}));
    addOptional(p, 'tractType', 'dwi', @(x) ismember(x, {'dwi', 'template'}));
    addParameter(p, 'force', false, @islogical);
    addParameter(p, 'test', false, @islogical);
    addParameter(p, 'lesionMask', 'manual', @(x) ismember(x, {'LSTAI', 'manual'}));
    parse(p, varargin{:});

    % Load configuration
    config = load_config();

    try
        % Setup test environment if needed
        if p.Results.test
            config = setupTestEnvironment(config, ['gen_' p.Results.metric '_metrics']);
        end

        % Initialize logging
        logHelper(config, 'initialize');

        % Setup output file path using config
        if isfield(config, 'isTest') && config.isTest
            % Modify output filename to include "All"
            origName = config.excelOutputs.(p.Results.tractType).(lower(p.Results.metric));
            [~, name, ext] = fileparts(origName);
            newName = [name '_All' ext];
            outputFile = fullfile(config.testDir, newName);
        else
            % Modify output filename to include "All"
            origName = config.excelOutputs.(p.Results.tractType).(lower(p.Results.metric));
            [~, name, ext] = fileparts(origName);
            newName = [name '_All' ext];
            outputFile = fullfile(config.bidsDir, config.statsDir, newName);
        end

        % Create output directory if it doesn't exist
        outputDir = fileparts(outputFile);
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end

        % Initialize or load existing data table
        dataTable = initializeDataTable(outputFile);

        % Get subject list (including controls)
        subjectIDs = getSubjectList(config, true);  % true to include controls

        % Define tract groups
        groups = {'Association', 'Cerebellar', 'Occipitoparietal', 'ProjectionBrainstem'};

        % Process subjects
        processAllSubjects(subjectIDs, dataTable, config, groups, outputFile, p.Results.force, p.Results.metric, p.Results.tractType, p.Results.lesionMask);

    catch ME
        logError(ME);
        rethrow(ME);
    end

    % Finalize logging
    logHelper([], 'finalize');
end

function subjectIDs = getSubjectList(config, includeControls)
    % Get all subject directories
    allSubjects = dir(fullfile(config.bidsDir, config.derivDir, 'TractoFlow_post', 'sub-*'));
    allSubjects = allSubjects([allSubjects.isdir]);

    if isempty(allSubjects)
        allSubjects = dir(fullfile(config.bidsDir, config.derivDir, 'maps', 'sub-*'));
        allSubjects = allSubjects([allSubjects.isdir]);
    end

    % Filter controls if needed
    if ~includeControls
        isControl = startsWith({allSubjects.name}, 'sub-C');
        allSubjects = allSubjects(~isControl);
    end

    % In test mode, limit to first few subjects
    if isfield(config, 'isTest') && config.isTest
        allSubjects = allSubjects(1:min(3, length(allSubjects)));
    end

    subjectIDs = {allSubjects.name};

    if isempty(subjectIDs)
        error('No subjects found to process');
    end
end

function processAllSubjects(subjectIDs, dataTable, config, groups, outputFile, force, metric, tractType, lesionMask)
    % Initialize results container
    T_cell = cell(length(subjectIDs), 1);

    % Define all expected columns (now with _All, and _Tail in addition to _Lesion and _NAWM)
    expectedColumns = {'SubjectID'};
    for iGroup = 1:length(groups)
        group = groups{iGroup};
        expectedColumns = [expectedColumns, ...
            {[group 'L_Lesion'], [group 'L_NAWM'], [group 'L_All'], [group 'L_Tail'], ...
             [group 'R_Lesion'], [group 'R_NAWM'], [group 'R_All'], [group 'R_Tail'], ...
             [group '_Lesion'], [group '_NAWM'], [group '_All'], [group '_Tail']}];
    end

    % Start parallel pool if not already running
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        fprintf('Starting parallel pool with %d workers...\n', maxNumCompThreads);
        parpool('local', maxNumCompThreads);
    else
        fprintf('Using existing parallel pool with %d workers\n', poolobj.NumWorkers);
    end

    % Process subjects in parallel
    fprintf('Processing %d subjects...\n', length(subjectIDs));
    tic;

    parfor i = 1:length(subjectIDs)
        subjectID = subjectIDs{i};
        T = processSingleSubject(subjectID, config, groups, metric, tractType, lesionMask);

        % Create a new table with all expected columns
        newT = table('Size', [1, length(expectedColumns)], ...
                    'VariableTypes', [{'cell'}, repmat({'double'}, 1, length(expectedColumns)-1)], ...
                    'VariableNames', expectedColumns);

        % Set subject ID
        newT.SubjectID = {subjectID};

        % Copy existing values and fill missing ones with NaN
        existingCols = T.Properties.VariableNames;
        for j = 1:length(expectedColumns)
            colName = expectedColumns{j};
            if ismember(colName, existingCols)
                newT.(colName) = T.(colName);
            else
                newT.(colName) = NaN;
            end
        end

        T_cell{i} = newT;
    end

    elapsed = toc;
    fprintf('Completed processing %d subjects in %.1f seconds (%.2f sec/subject)\n', ...
        length(subjectIDs), elapsed, elapsed/length(subjectIDs));

    % Combine results and save
    newDataTable = vertcat(T_cell{:});
    if ~isempty(dataTable)
        dataTable = outerjoin(dataTable, newDataTable, 'Keys', 'SubjectID', 'MergeKeys', true);
    else
        dataTable = newDataTable;
    end

    % Save results
    writetable(dataTable, outputFile);
    fprintf('Results saved to: %s\n', outputFile);
end

function T = processSingleSubject(subjectID, config, groups, metric, tractType, lesionMask)
    % fprintf('Processing %s...\n', subjectID);  % Commented out to reduce verbosity
    isControl = startsWith(subjectID, 'sub-C');

    % Initialize table with subject ID
    T = table({subjectID}, 'VariableNames', {'SubjectID'});

    % Load metric map
    switch metric
        case 'T1'
            mapPath = fullfile(config.bidsDir, config.derivDir, 'maps', subjectID, ...
                             [subjectID '_ses-01_space-individual_T1map.nii.gz']);
        case 'FA'
            mapPath = fullfile(config.bidsDir, config.derivDir, 'maps', subjectID, ...
                             [subjectID '_space-individual_FA.nii.gz']);
        case 'MD'
            mapPath = fullfile(config.bidsDir, config.derivDir, 'maps', subjectID, ...
                             [subjectID '_space-individual_MD.nii.gz']);
        case 'MTR'
            mapPath = fullfile(config.bidsDir, config.derivDir, 'maps', subjectID, ...
                             [subjectID '_ses-01_space-individual_MTRmap.nii.gz']);
        case 'AD'
            mapPath = fullfile(config.bidsDir, config.derivDir, 'maps', subjectID, ...
                             [subjectID '_space-individual_AD.nii.gz']);
        case 'RD'
            mapPath = fullfile(config.bidsDir, config.derivDir, 'maps', subjectID, ...
                             [subjectID '_space-individual_RD.nii.gz']);
    end

    if ~exist(mapPath, 'file')
        warning('%s map not found for %s', metric, subjectID);
        return;
    end
    metricMap = niftiread(mapPath);

    % Load lesion mask for MS patients
    lesionMaskData = loadLesionMask(config, subjectID, lesionMask);

    % Process each tract group
    for iGroup = 1:length(groups)
        group = groups{iGroup};

        % Process left hemisphere
        [lesionValL, nawmValL, allValL, tailValL] = processHemisphereTract(config, subjectID, group, 'L', ...
                                                       metricMap, lesionMaskData, isControl, tractType, metric);
        T.([group 'L_Lesion']) = lesionValL;
        T.([group 'L_NAWM']) = nawmValL;
        T.([group 'L_All']) = allValL;
        T.([group 'L_Tail']) = tailValL;

        % Process right hemisphere
        [lesionValR, nawmValR, allValR, tailValR] = processHemisphereTract(config, subjectID, group, 'R', ...
                                                       metricMap, lesionMaskData, isControl, tractType, metric);
        T.([group 'R_Lesion']) = lesionValR;
        T.([group 'R_NAWM']) = nawmValR;
        T.([group 'R_All']) = allValR;
        T.([group 'R_Tail']) = tailValR;

         % Process whole brain (both hemispheres combined)
        [lesionVal, nawmVal, allVal, tailVal] = processHemisphereTract(config, subjectID, group, '', ...
                    metricMap, lesionMaskData, isControl, tractType, metric);
        T.([group '_Lesion']) = lesionVal;
        T.([group '_NAWM']) = nawmVal;
        T.([group '_All']) = allVal;
        T.([group '_Tail']) = tailVal;
    end
end

function [lesionVal, nawmVal, allVal, tailVal] = processHemisphereTract(config, subjectID, group, hemi, ...
                                                      metricMap, lesionMask, isControl, tractType, metric)
    % Initialize outputs
    lesionVal = 0;
    nawmVal = 0;
    allVal = 0;
    tailVal = 0;

    % Load tract mask based on tract type
    if strcmp(tractType, 'template')
        tractPath = fullfile(config.bidsDir, config.derivDir, 'TractoFlow_post', ...
                            subjectID, 'groupTractTemplate', [group hemi '.nii.gz']);
    else
        tractPath = fullfile(config.bidsDir, config.derivDir, 'TractoFlow_post', ...
                            subjectID, 'groupTract', [group hemi '.nii.gz']);
    end

    if ~exist(tractPath, 'file')
        warning('Tract mask not found: %s', tractPath);
        return;
    end

    % Load and reorient tract mask (x -z -y)
    tractMask = niftiread(tractPath);
    tractMask = flip(flip(permute(tractMask, [1 3 2]), 2), 3);  % Combined permute and flips
    tractMask = logical(tractMask);

    % Early exit if no tract voxels
    if ~any(tractMask(:))
        return;
    end

    % Extract all values from tract once (more efficient than multiple indexing)
    allValues = metricMap(tractMask);

    % Remove NaN/invalid values upfront
    validIdx = isfinite(allValues);
    allValues_clean = allValues(validIdx);

    if isempty(allValues_clean)
        return;
    end

    if isControl
        % For controls: no lesions, so NAWM = All = entire tract
        nawmVal = mean(allValues_clean);
        allVal = nawmVal;
        tailVal = calculateTailPercentile(allValues_clean, metric);
    else
        % For MS patients, separate lesions, NAWM, and All
        if ~isempty(lesionMask)
            % Calculate All first (we already have the values)
            allVal = mean(allValues_clean);
            tailVal = calculateTailPercentile(allValues_clean, metric);

            % Get lesion values
            lesionMaskInTract = tractMask & lesionMask;
            if any(lesionMaskInTract(:))
                lesionValues = metricMap(lesionMaskInTract);
                lesionValues = lesionValues(isfinite(lesionValues));
                if ~isempty(lesionValues)
                    lesionVal = mean(lesionValues);
                end
            end

            % Get NAWM values
            nawmMask = tractMask & ~lesionMask;
            if any(nawmMask(:))
                nawmValues = metricMap(nawmMask);
                nawmValues = nawmValues(isfinite(nawmValues));
                if ~isempty(nawmValues)
                    nawmVal = mean(nawmValues);
                end
            end
        else
            % No lesion mask available for MS patient - treat as all NAWM
            nawmVal = mean(allValues_clean);
            allVal = nawmVal;
            tailVal = calculateTailPercentile(allValues_clean, metric);
        end
    end
end

function tailVal = calculateTailPercentile(values, metric)
    % Calculate tail percentile based on metric type
    % MTR, FA: use 10th percentile (lower tail = worse)
    % T1, MD, AD, RD: use 90th percentile (upper tail = worse)
    %
    % NOTE: Input values should already have NaN/Inf removed for efficiency

    if isempty(values)
        tailVal = 0;
        return;
    end

    % Determine which percentile to use based on metric
    switch metric
        case {'MTR', 'FA'}
            % Lower values indicate damage - use 10th percentile
            tailVal = prctile(values, 10);
        case {'T1', 'MD', 'AD', 'RD'}
            % Higher values indicate damage - use 90th percentile
            tailVal = prctile(values, 90);
        otherwise
            % Default to median
            tailVal = median(values);
    end
end

function dataTable = initializeDataTable(outputFile)
    if exist(outputFile, 'file')
        response = input(['Output file ' outputFile ' already exists. Overwrite? (y/n): '], 's');
        if strcmpi(response, 'y') || strcmpi(response, 'yes')
            dataTable = table();
        else
            error('Operation cancelled by user.');
        end
    else
        dataTable = table();
    end
end

function [lesionMask] = loadLesionMask(config, subjectID, lesionMaskType)
    % Construct lesion mask path based on type
    if strcmp(lesionMaskType, 'LSTAI')
        lesionSegPath = fullfile(config.bidsDir, config.derivDir, 'lesionMask', ...
                               subjectID, 'ses-01', [subjectID '_ses-01_desc-LSTAI_mask.nii.gz']);
    else % manual
        lesionSegPath = fullfile(config.bidsDir, config.derivDir, 'lesionMask', ...
                               subjectID, 'ses-01', [subjectID '_ses-01_desc-lesionManual_mask.nii.gz']);
    end

    if exist(lesionSegPath, 'file')
        try
            lesionMask = double(niftiread(lesionSegPath));
        catch ME
            warning(ME.identifier, '%s', ME.message);
            lesionMask = zeros(256, 256, 256);  % Default size
        end
    else
        warning('Lesion segmentation file not found for %s: %s', subjectID, lesionSegPath);
        lesionMask = zeros(256, 256, 256);  % Default size
    end
end
