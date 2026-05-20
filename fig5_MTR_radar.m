function fig5_MTR_radar(varargin)
% F9 Creates polar plots of MTR values in lesions and NAWM by MS type
%
% Usage:
%   f9()          % Uses default DWI group tracts, both hemispheres
%   f9('dwi')     % Uses DWI group tracts, both hemispheres
%   f9('template') % Uses template group tracts, both hemispheres
%   f9('dwi', 'hemi')     % Uses DWI group tracts, separate hemispheres
%   f9('template', 'hemi') % Uses template group tracts, separate hemispheres
%
% This script creates polar plots:
% - With 'hemi' option: Four polar plots (2x2 layout)
%   1. Top left: Lesion MTR values (left hemisphere)
%   2. Top right: Lesion MTR values (right hemisphere)
%   3. Bottom left: NAWM MTR values (left hemisphere)
%   4. Bottom right: NAWM MTR values (right hemisphere)
% - Without 'hemi' option: Two polar plots (2x1 layout)
%   1. Top: Lesion MTR values (combined hemispheres)
%   2. Bottom: NAWM MTR values (combined hemispheres)
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
                        config.excelOutputs.(tractType).mtr);
    
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
    
    % Initialize arrays for MTR values
    numGroups = length(groups);
    if separateHemi
        lesionMTRL = zeros(numGroups, 3); % 3 MS types (excluding controls)
        lesionMTRR = zeros(numGroups, 3);
        nawmMTRL = zeros(numGroups, 4);   % 4 types (including controls)
        nawmMTRR = zeros(numGroups, 4);
        
        % Initialize arrays for standard errors
        lesionMTRL_se = zeros(numGroups, 3);
        lesionMTRR_se = zeros(numGroups, 3);
        nawmMTRL_se = zeros(numGroups, 4);
        nawmMTRR_se = zeros(numGroups, 4);
    else
        lesionMTR = zeros(numGroups, 3);
        nawmMTR = zeros(numGroups, 4);
        lesionMTR_se = zeros(numGroups, 3);
        nawmMTR_se = zeros(numGroups, 4);
    end
    
    % Calculate averages for each group and type
    msTypesList = {'RRMS', 'SPMS', 'PPMS'};
    allTypesList = [msTypesList, {'control'}];
    
    for i = 1:numGroups
        group = groups{i};
        
        % Process lesion MTR (MS patients only)
        for j = 1:length(msTypesList)
            msType = msTypesList{j};
            idx = data.MSType == msType;
            
            if separateHemi
                % Left hemisphere lesions
                colL = [group 'L_Lesion'];
                if ismember(colL, data.Properties.VariableNames)
                    values = data{idx, colL};
                    values = values(values ~= 0);
                    lesionMTRL(i,j) = mean(values, 'omitnan');
                    lesionMTRL_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
                
                % Right hemisphere lesions
                colR = [group 'R_Lesion'];
                if ismember(colR, data.Properties.VariableNames)
                    values = data{idx, colR};
                    values = values(values ~= 0);
                    lesionMTRR(i,j) = mean(values, 'omitnan');
                    lesionMTRR_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            else
                % Combined hemisphere lesions
                col = [group '_Lesion'];
                if ismember(col, data.Properties.VariableNames)
                    values = data{idx, col};
                    values = values(values ~= 0);
                    lesionMTR(i,j) = mean(values, 'omitnan');
                    lesionMTR_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            end
        end
        
        % Process NAWM MTR (all subjects)
        for j = 1:length(allTypesList)
            type = allTypesList{j};
            idx = data.MSType == type;
            
            if separateHemi
                % Left hemisphere NAWM
                colL = [group 'L_NAWM'];
                if ismember(colL, data.Properties.VariableNames)
                    values = data{idx, colL};
                    values = values(values ~= 0);
                    nawmMTRL(i,j) = mean(values, 'omitnan');
                    nawmMTRL_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
                
                % Right hemisphere NAWM
                colR = [group 'R_NAWM'];
                if ismember(colR, data.Properties.VariableNames)
                    values = data{idx, colR};
                    values = values(values ~= 0);
                    nawmMTRR(i,j) = mean(values, 'omitnan');
                    nawmMTRR_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            else
                % Combined hemisphere NAWM
                col = [group '_NAWM'];
                if ismember(col, data.Properties.VariableNames)
                    values = data{idx, col};
                    values = values(values ~= 0);
                    nawmMTR(i,j) = mean(values, 'omitnan');
                    nawmMTR_se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            end
        end
    end
    
    % Create the figure
    if separateHemi
        fig = figure('Position', [100 100 1200 1000]);
        maxCommonR = max(max([lesionMTRL(:) + lesionMTRL_se(:); 
                             lesionMTRR(:) + lesionMTRR_se(:);
                             nawmMTRL(:) + nawmMTRL_se(:);
                             nawmMTRR(:) + nawmMTRR_se(:)])) * 1.1;
    else
        fig = figure('Position', [100 100 600 1000]);
        maxCommonR = max(max([lesionMTR(:) + lesionMTR_se(:);
                             nawmMTR(:) + nawmMTR_se(:)])) * 1.1;
    end
    
    if maxCommonR == 0
        maxCommonR = 1;
    end
    
    if separateHemi
        % Plot lesion MTR values (left hemisphere)
        subplot(2,2,1);
        createPolarPlot(lesionMTRL, lesionMTRL_se, groups, msTypesList, config, ...
                        'Lesion MTR (Left)', tractType, maxCommonR, false);
        
        % Plot lesion MTR values (right hemisphere)
        subplot(2,2,2);
        createPolarPlot(lesionMTRR, lesionMTRR_se, groups, msTypesList, config, ...
                        'Lesion MTR (Right)', tractType, maxCommonR, false);
        
        % Plot NAWM MTR values (left hemisphere)
        subplot(2,2,3);
        createPolarPlot(nawmMTRL, nawmMTRL_se, groups, allTypesList, config, ...
                        'NAWM MTR (Left)', tractType, maxCommonR, true);
        
        % Plot NAWM MTR values (right hemisphere)
        subplot(2,2,4);
        createPolarPlot(nawmMTRR, nawmMTRR_se, groups, allTypesList, config, ...
                        'NAWM MTR (Right)', tractType, maxCommonR, true);
    else
        % Plot lesion MTR values (combined hemispheres)
        subplot(2,1,1);
        createPolarPlot(lesionMTR, lesionMTR_se, groups, msTypesList, config, ...
                        'Lesion MTR', tractType, maxCommonR, false);
        
        % Plot NAWM MTR values (combined hemispheres)
        subplot(2,1,2);
        createPolarPlot(nawmMTR, nawmMTR_se, groups, allTypesList, config, ...
                        'NAWM MTR', tractType, maxCommonR, true);
    end
    
    % Add a common title
    if strcmp(tractType, 'template')
        titleStr = 'MTR Values in Lesions and NAWM by MS Type (Atlas)';
    else
        titleStr = 'MTR Values in Lesions and NAWM by MS Type';
    end
    sgtitle(titleStr, 'FontSize', 14);
end
