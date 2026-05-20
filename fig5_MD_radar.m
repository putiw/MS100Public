function fig5_MD_radar(varargin)
% F6 Creates polar plots of MD values in lesions and NAWM by MS type
%
% Usage:
%   f6()          % Uses default DWI group tracts, both hemispheres
%   f6('dwi')     % Uses DWI group tracts, both hemispheres
%   f6('template') % Uses template group tracts, both hemispheres
%   f6('dwi', 'hemi')     % Uses DWI group tracts, separate hemispheres
%   f6('template', 'hemi') % Uses template group tracts, separate hemispheres
%
% This script creates polar plots:
% - With 'hemi' option: Four polar plots (2x2 layout)
%   1. Top left: Lesion MD values (left hemisphere)
%   2. Top right: Lesion MD values (right hemisphere)
%   3. Bottom left: NAWM MD values (left hemisphere)
%   4. Bottom right: NAWM MD values (right hemisphere)
% - Without 'hemi' option: Two polar plots (2x1 layout)
%   1. Top: Lesion MD values (combined hemispheres)
%   2. Bottom: NAWM MD values (combined hemispheres)
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
                        config.excelOutputs.(tractType).md);
    
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
    
    % Initialize arrays for MD values
    numGroups = length(groups);
    if separateHemi
        lesionMDL = zeros(numGroups, 3); % 3 MS types (excluding controls)
        lesionMDR = zeros(numGroups, 3);
        nawmMDL = zeros(numGroups, 4);   % 4 types (including controls)
        nawmMDR = zeros(numGroups, 4);
        
        % Initialize arrays for standard errors
        lesionMDL_se = zeros(numGroups, 3);
        lesionMDR_se = zeros(numGroups, 3);
        nawmMDL_se = zeros(numGroups, 4);
        nawmMDR_se = zeros(numGroups, 4);
    else
        lesionMD = zeros(numGroups, 3);
        nawmMD = zeros(numGroups, 4);
        lesionMD_se = zeros(numGroups, 3);
        nawmMD_se = zeros(numGroups, 4);
    end
    
    % Calculate averages for each group and type
    msTypesList = {'RRMS', 'SPMS', 'PPMS'};
    allTypesList = [msTypesList, {'control'}];
    
    for i = 1:numGroups
        group = groups{i};
        
        % Process lesion MD (MS patients only)
        for j = 1:length(msTypesList)
            msType = msTypesList{j};
            idx = data.MSType == msType;
            
            if separateHemi
                % Left hemisphere lesions
                colL = [group 'L_Lesion'];
                if ismember(colL, data.Properties.VariableNames)
                    values = data{idx, colL};
                    values = values(values ~= 0);
                    lesionMDL(i,j) = mean(values, 'omitnan');
                    lesionMDL_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
                
                % Right hemisphere lesions
                colR = [group 'R_Lesion'];
                if ismember(colR, data.Properties.VariableNames)
                    values = data{idx, colR};
                    values = values(values ~= 0);
                    lesionMDR(i,j) = mean(values, 'omitnan');
                    lesionMDR_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            else
                % Combined hemisphere lesions
                col = [group '_Lesion'];
                if ismember(col, data.Properties.VariableNames)
                    values = data{idx, col};
                    values = values(values ~= 0);
                    lesionMD(i,j) = mean(values, 'omitnan');
                    lesionMD_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            end
        end
        
        % Process NAWM MD (all subjects)
        for j = 1:length(allTypesList)
            type = allTypesList{j};
            idx = data.MSType == type;
            
            if separateHemi
                % Left hemisphere NAWM
                colL = [group 'L_NAWM'];
                if ismember(colL, data.Properties.VariableNames)
                    values = data{idx, colL};
                    values = values(values ~= 0);
                    nawmMDL(i,j) = mean(values, 'omitnan');
                    nawmMDL_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
                
                % Right hemisphere NAWM
                colR = [group 'R_NAWM'];
                if ismember(colR, data.Properties.VariableNames)
                    values = data{idx, colR};
                    values = values(values ~= 0);
                    nawmMDR(i,j) = mean(values, 'omitnan');
                    nawmMDR_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            else
                % Combined hemisphere NAWM
                col = [group '_NAWM'];
                if ismember(col, data.Properties.VariableNames)
                    values = data{idx, col};
                    values = values(values ~= 0);
                    nawmMD(i,j) = mean(values, 'omitnan');
                    nawmMD_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            end
        end
    end
    
    % Create the figure
    if separateHemi
        fig = figure('Position', [100 100 1200 1000]);
        maxCommonR = max(max([lesionMDL(:) + lesionMDL_se(:); 
                             lesionMDR(:) + lesionMDR_se(:);
                             nawmMDL(:) + nawmMDL_se(:);
                             nawmMDR(:) + nawmMDR_se(:)])) * 1.1;
    else
        fig = figure('Position', [100 100 600 1000]);
        maxCommonR = max(max([lesionMD(:) + lesionMD_se(:);
                             nawmMD(:) + nawmMD_se(:)])) * 1.1;
    end
    
    if maxCommonR == 0
        maxCommonR = 1;
    end
    
    if separateHemi
        % Plot lesion MD values (left hemisphere)
        subplot(2,2,1);
        createPolarPlot(lesionMDL, lesionMDL_se, groups, msTypesList, config, ...
                        'Lesion MD (Left)', tractType, maxCommonR, false);
        
        % Plot lesion MD values (right hemisphere)
        subplot(2,2,2);
        createPolarPlot(lesionMDR, lesionMDR_se, groups, msTypesList, config, ...
                        'Lesion MD (Right)', tractType, maxCommonR, false);
        
        % Plot NAWM MD values (left hemisphere)
        subplot(2,2,3);
        createPolarPlot(nawmMDL, nawmMDL_se, groups, allTypesList, config, ...
                        'NAWM MD (Left)', tractType, maxCommonR, true);
        
        % Plot NAWM MD values (right hemisphere)
        subplot(2,2,4);
        createPolarPlot(nawmMDR, nawmMDR_se, groups, allTypesList, config, ...
                        'NAWM MD (Right)', tractType, maxCommonR, true);
    else
        % Plot lesion MD values (combined hemispheres)
        subplot(2,1,1);
        createPolarPlot(lesionMD, lesionMD_se, groups, msTypesList, config, ...
                        'Lesion MD', tractType, maxCommonR, false);
        
        % Plot NAWM MD values (combined hemispheres)
        subplot(2,1,2);
        createPolarPlot(nawmMD, nawmMD_se, groups, allTypesList, config, ...
                        'NAWM MD', tractType, maxCommonR, true);
    end
    
    % Add a common title
    if strcmp(tractType, 'template')
        titleStr = 'MD Values in Lesions and NAWM by MS Type (Atlas)';
    else
        titleStr = 'MD Values in Lesions and NAWM by MS Type';
    end
    sgtitle(titleStr, 'FontSize', 14);
end
