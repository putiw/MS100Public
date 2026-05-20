function fig5_T1_radar(varargin)
% F4 Creates polar plots of T1 values in lesions and NAWM by MS type
%
% Usage:
%   f4()          % Uses default DWI group tracts, both hemispheres
%   f4('dwi')     % Uses DWI group tracts, both hemispheres
%   f4('template') % Uses template group tracts, both hemispheres
%   f4('dwi', 'hemi')     % Uses DWI group tracts, separate hemispheres
%   f4('template', 'hemi') % Uses template group tracts, separate hemispheres
%
% This script creates polar plots:
% - With 'hemi' option: Four polar plots (2x2 layout)
%   1. Top left: Lesion T1 values (left hemisphere)
%   2. Top right: Lesion T1 values (right hemisphere)
%   3. Bottom left: NAWM T1 values (left hemisphere)
%   4. Bottom right: NAWM T1 values (right hemisphere)
% - Without 'hemi' option: Two polar plots (2x1 layout)
%   1. Top: Lesion T1 values (combined hemispheres)
%   2. Bottom: NAWM T1 values (combined hemispheres)
%
% Groups plotted:
% - Association
% - Cerebellar
% - Occipitoparietal
% - ProjectionBrainstem

    % Parse input arguments
    p = inputParser;
    addOptional(p, 'tractType', 'dwi', @(x) ismember(x, {'dwi', 'template'}));
    addOptional(p, 'hemi', '', @(x) strcmp(x, 'hemi'));
    parse(p, varargin{:});
    tractType = p.Results.tractType;
    separateHemi = strcmp(p.Results.hemi, 'hemi');

    % Add helper functions to path
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));
    
    % Load configuration
    config = load_config();
    
    % Get the correct Excel file path from config
    excelFile = fullfile(config.bidsDir, config.statsDir, ...
                        config.excelOutputs.(tractType).t1);
    
    if ~exist(excelFile, 'file')
        error('Excel file not found: %s', excelFile);
    end
    data = readtable(excelFile);
    
    % Extract MS types from subject IDs
    msTypes = cell(height(data), 1);
    for i = 1:height(data)
        subID = data.SubjectID{i};
        if contains(subID, 'PP')
            msTypes{i} = 'PPMS';
        elseif contains(subID, 'RR')
            msTypes{i} = 'RRMS';
        elseif contains(subID, 'SP')
            msTypes{i} = 'SPMS';
        elseif contains(subID, 'C')
            msTypes{i} = 'control';
        else
            warning('Unknown subject type for subject: %s', subID);
            msTypes{i} = 'Unknown';
        end
    end
    data.MSType = categorical(msTypes);
    
    % Define the groups
    groups = {'Association', 'Cerebellar', 'Occipitoparietal', 'ProjectionBrainstem'};
    
    % Initialize arrays for T1 values
    numGroups = length(groups);
    if separateHemi
        lesionT1L = zeros(numGroups, 3); % 3 MS types (excluding controls)
        lesionT1R = zeros(numGroups, 3);
        nawmT1L = zeros(numGroups, 4);   % 4 types (including controls)
        nawmT1R = zeros(numGroups, 4);
        
        % Initialize arrays for standard errors
        lesionT1L_se = zeros(numGroups, 3);
        lesionT1R_se = zeros(numGroups, 3);
        nawmT1L_se = zeros(numGroups, 4);
        nawmT1R_se = zeros(numGroups, 4);
    else
        lesionT1 = zeros(numGroups, 3);
        nawmT1 = zeros(numGroups, 4);
        lesionT1_se = zeros(numGroups, 3);
        nawmT1_se = zeros(numGroups, 4);
    end
    
    % Calculate averages for each group and type
    msTypesList = {'RRMS', 'SPMS', 'PPMS'};
    allTypesList = [msTypesList, {'control'}];
    
    % Helper function to process values and calculate statistics
    function [meanVal, seVal] = processValues(data, idx, colName)
        if ismember(colName, data.Properties.VariableNames)
            values = data{idx, colName};
            values = values(values ~= 0);
            meanVal = mean(values, 'omitnan');
            seVal = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
        else
            meanVal = NaN;
            seVal = NaN;
        end
    end

    % Helper function to get column name based on group and type
    function colName = getColumnName(group, type, hemi)
        if nargin < 3
            hemi = '';
        end
        colName = [group hemi '_' type];
    end

    % Process values for each group and type
    for i = 1:numGroups
        group = groups{i};
        
        % Process lesion T1 (MS patients only)
        for j = 1:length(msTypesList)
            msType = msTypesList{j};
            idx = data.MSType == msType;
            
            if separateHemi
                % Left hemisphere lesions
                [lesionT1L(i,j), lesionT1L_se(i,j)] = processValues(data, idx, getColumnName(group, 'Lesion', 'L'));
                
                % Right hemisphere lesions
                [lesionT1R(i,j), lesionT1R_se(i,j)] = processValues(data, idx, getColumnName(group, 'Lesion', 'R'));
            else
                % Combined hemisphere lesions
                [lesionT1(i,j), lesionT1_se(i,j)] = processValues(data, idx, getColumnName(group, 'Lesion'));
            end
        end
        
        % Process NAWM T1 (all subjects)
        for j = 1:length(allTypesList)
            type = allTypesList{j};
            idx = data.MSType == type;
            
            if separateHemi
                % Left hemisphere NAWM
                [nawmT1L(i,j), nawmT1L_se(i,j)] = processValues(data, idx, getColumnName(group, 'NAWM', 'L'));
                
                % Right hemisphere NAWM
                [nawmT1R(i,j), nawmT1R_se(i,j)] = processValues(data, idx, getColumnName(group, 'NAWM', 'R'));
            else
                % Combined hemisphere NAWM
                [nawmT1(i,j), nawmT1_se(i,j)] = processValues(data, idx, getColumnName(group, 'NAWM'));
            end
        end
    end
    
    % Create the figure
    if separateHemi
        fig = figure('Position', [100 100 1200 1000]);
        maxCommonR = max(max([lesionT1L(:) + lesionT1L_se(:); 
                             lesionT1R(:) + lesionT1R_se(:);
                             nawmT1L(:) + nawmT1L_se(:);
                             nawmT1R(:) + nawmT1R_se(:)])) * 1.1;
    else
        fig = figure('Position', [100 100 600 1000]);
        maxCommonR = max(max([lesionT1(:) + lesionT1_se(:);
                             nawmT1(:) + nawmT1_se(:)])) * 1.1;
    end
    
    if maxCommonR == 0
        maxCommonR = 1;
    end
    
    if separateHemi
        % Plot lesion T1 values (left hemisphere)
        subplot(2,2,1);
        createPolarPlot(lesionT1L, lesionT1L_se, groups, msTypesList, config, ...
                        'Lesion T1 (Left)', tractType, maxCommonR, false);
        
        % Plot lesion T1 values (right hemisphere)
        subplot(2,2,2);
        createPolarPlot(lesionT1R, lesionT1R_se, groups, msTypesList, config, ...
                        'Lesion T1 (Right)', tractType, maxCommonR, false);
        
        % Plot NAWM T1 values (left hemisphere)
        subplot(2,2,3);
        createPolarPlot(nawmT1L, nawmT1L_se, groups, allTypesList, config, ...
                        'NAWM T1 (Left)', tractType, maxCommonR, true);
        
        % Plot NAWM T1 values (right hemisphere)
        subplot(2,2,4);
        createPolarPlot(nawmT1R, nawmT1R_se, groups, allTypesList, config, ...
                        'NAWM T1 (Right)', tractType, maxCommonR, true);
    else
        % Plot lesion T1 values (combined hemispheres)
        subplot(2,1,1);
        createPolarPlot(lesionT1, lesionT1_se, groups, msTypesList, config, ...
                        'Lesion T1', tractType, maxCommonR, false);
        
        % Plot NAWM T1 values (combined hemispheres)
        subplot(2,1,2);
        createPolarPlot(nawmT1, nawmT1_se, groups, allTypesList, config, ...
                        'NAWM T1', tractType, maxCommonR, true);
    end
    
    % Add a common title
    if strcmp(tractType, 'template')
        titleStr = 'T1 Values in Lesions and NAWM by MS Type (Atlas)';
    else
        titleStr = 'T1 Values in Lesions and NAWM by MS Type';
    end
    sgtitle(titleStr, 'FontSize', 14);
end

