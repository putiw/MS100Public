function fig5_FA_radar(varargin)
% F5 Creates polar plots of FA values in lesions and NAWM by MS type
%
% Usage:
%   f5()          % Uses default DWI group tracts, both hemispheres
%   f5('dwi')     % Uses DWI group tracts, both hemispheres
%   f5('template') % Uses template group tracts, both hemispheres
%   f5('dwi', 'hemi')     % Uses DWI group tracts, separate hemispheres
%   f5('template', 'hemi') % Uses template group tracts, separate hemispheres
%
% This script creates polar plots:
% - With 'hemi' option: Four polar plots (2x2 layout)
%   1. Top left: Lesion FA values (left hemisphere)
%   2. Top right: Lesion FA values (right hemisphere)
%   3. Bottom left: NAWM FA values (left hemisphere)
%   4. Bottom right: NAWM FA values (right hemisphere)
% - Without 'hemi' option: Two polar plots (2x1 layout)
%   1. Top: Lesion FA values (combined hemispheres)
%   2. Bottom: NAWM FA values (combined hemispheres)
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
                        config.excelOutputs.(tractType).fa);
    
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
    
    % Initialize arrays for FA values
    numGroups = length(groups);
    if separateHemi
        lesionFAL = zeros(numGroups, 3); % 3 MS types (excluding controls)
        lesionFAR = zeros(numGroups, 3);
        nawmFAL = zeros(numGroups, 4);   % 4 types (including controls)
        nawmFAR = zeros(numGroups, 4);
        
        % Initialize arrays for standard errors
        lesionFAL_se = zeros(numGroups, 3);
        lesionFAR_se = zeros(numGroups, 3);
        nawmFAL_se = zeros(numGroups, 4);
        nawmFAR_se = zeros(numGroups, 4);
    else
        lesionFA = zeros(numGroups, 3);
        nawmFA = zeros(numGroups, 4);
        lesionFA_se = zeros(numGroups, 3);
        nawmFA_se = zeros(numGroups, 4);
    end
    
    % Calculate averages for each group and type
    msTypesList = {'RRMS', 'SPMS', 'PPMS'};
    allTypesList = [msTypesList, {'control'}];
    
    for i = 1:numGroups
        group = groups{i};
        
        % Process lesion FA (MS patients only)
        for j = 1:length(msTypesList)
            msType = msTypesList{j};
            idx = data.MSType == msType;
            
            if separateHemi
                % Left hemisphere lesions
                colL = [group 'L_Lesion'];
                if ismember(colL, data.Properties.VariableNames)
                    values = data{idx, colL};
                    values = values(values ~= 0);
                    lesionFAL(i,j) = mean(values, 'omitnan');
                    lesionFAL_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
                
                % Right hemisphere lesions
                colR = [group 'R_Lesion'];
                if ismember(colR, data.Properties.VariableNames)
                    values = data{idx, colR};
                    values = values(values ~= 0);
                    lesionFAR(i,j) = mean(values, 'omitnan');
                    lesionFAR_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            else
                % Combined hemisphere lesions
                col = [group '_Lesion'];
                if ismember(col, data.Properties.VariableNames)
                    values = data{idx, col};
                    values = values(values ~= 0);
                    lesionFA(i,j) = mean(values, 'omitnan');
                    lesionFA_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            end
        end
        
        % Process NAWM FA (all subjects)
        for j = 1:length(allTypesList)
            type = allTypesList{j};
            idx = data.MSType == type;
            
            if separateHemi
                % Left hemisphere NAWM
                colL = [group 'L_NAWM'];
                if ismember(colL, data.Properties.VariableNames)
                    values = data{idx, colL};
                    values = values(values ~= 0);
                    nawmFAL(i,j) = mean(values, 'omitnan');
                    nawmFAL_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
                
                % Right hemisphere NAWM
                colR = [group 'R_NAWM'];
                if ismember(colR, data.Properties.VariableNames)
                    values = data{idx, colR};
                    values = values(values ~= 0);
                    nawmFAR(i,j) = mean(values, 'omitnan');
                    nawmFAR_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            else
                % Combined hemisphere NAWM
                col = [group '_NAWM'];
                if ismember(col, data.Properties.VariableNames)
                    values = data{idx, col};
                    values = values(values ~= 0);
                    nawmFA(i,j) = mean(values, 'omitnan');
                    nawmFA_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            end
        end
    end
    
    % Create the figure
    if separateHemi
        fig = figure('Position', [100 100 1200 1000]);
        maxCommonR = max(max([lesionFAL(:) + lesionFAL_se(:); 
                             lesionFAR(:) + lesionFAR_se(:);
                             nawmFAL(:) + nawmFAL_se(:);
                             nawmFAR(:) + nawmFAR_se(:)])) * 1.1;
    else
        fig = figure('Position', [100 100 600 1000]);
        maxCommonR = max(max([lesionFA(:) + lesionFA_se(:);
                             nawmFA(:) + nawmFA_se(:)])) * 1.1;
    end
    
    if maxCommonR == 0
        maxCommonR = 1;
    end
    
    if separateHemi
        % Plot lesion FA values (left hemisphere)
        subplot(2,2,1);
        createPolarPlot(lesionFAL, lesionFAL_se, groups, msTypesList, config, ...
                        'Lesion FA (Left)', tractType, maxCommonR, false);
        
        % Plot lesion FA values (right hemisphere)
        subplot(2,2,2);
        createPolarPlot(lesionFAR, lesionFAR_se, groups, msTypesList, config, ...
                        'Lesion FA (Right)', tractType, maxCommonR, false);
        
        % Plot NAWM FA values (left hemisphere)
        subplot(2,2,3);
        createPolarPlot(nawmFAL, nawmFAL_se, groups, allTypesList, config, ...
                        'NAWM FA (Left)', tractType, maxCommonR, true);
        
        % Plot NAWM FA values (right hemisphere)
        subplot(2,2,4);
        createPolarPlot(nawmFAR, nawmFAR_se, groups, allTypesList, config, ...
                        'NAWM FA (Right)', tractType, maxCommonR, true);
    else
        % Plot lesion FA values (combined hemispheres)
        subplot(2,1,1);
        createPolarPlot(lesionFA, lesionFA_se, groups, msTypesList, config, ...
                        'Lesion FA', tractType, maxCommonR, false);
        
        % Plot NAWM FA values (combined hemispheres)
        subplot(2,1,2);
        createPolarPlot(nawmFA, nawmFA_se, groups, allTypesList, config, ...
                        'NAWM FA', tractType, maxCommonR, true);
    end
    
    % Add a common title
    if strcmp(tractType, 'template')
        titleStr = 'FA Values in Lesions and NAWM by MS Type (Atlas)';
    else
        titleStr = 'FA Values in Lesions and NAWM by MS Type';
    end
    sgtitle(titleStr, 'FontSize', 14);
end

