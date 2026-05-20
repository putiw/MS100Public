function fig4_lesionload_tract_radar(varargin)
    % F3 Creates polar plots of lesion metrics by MS type for left and right hemispheres
    %
    % Usage:
    %   f3()                    % Uses default DWI group tracts, plots lesion numbers, both hemispheres
    %   f3('dwi')              % Uses DWI group tracts, plots lesion numbers, both hemispheres
    %   f3('template')         % Uses template group tracts, plots lesion numbers, both hemispheres
    %   f3('dwi', 'number')    % Uses DWI group tracts, plots lesion numbers, both hemispheres
    %   f3('dwi', 'volume')    % Uses DWI group tracts, plots lesion volumes, both hemispheres
    %   f3('dwi', 'norm')      % Uses DWI group tracts, plots normalized volumes, both hemispheres
    %   f3('dwi', 'number', 'hemi') % Uses DWI group tracts, plots lesion numbers, separate hemispheres
    %   f3('template', 'volume', 'hemi') % Uses template group tracts, plots lesion volumes, separate hemispheres
    %   f3('template', 'norm', 'hemi')   % Uses template group tracts, plots normalized volumes, separate hemispheres
    %   f3('number')           % Uses default DWI group tracts, plots lesion numbers, both hemispheres
    %   f3('volume')           % Uses default DWI group tracts, plots lesion volumes, both hemispheres
    %   f3('norm')             % Uses default DWI group tracts, plots normalized volumes, both hemispheres
    %   f3('number', 'hemi')   % Uses default DWI group tracts, plots lesion numbers, separate hemispheres
    %   f3('volume', 'hemi')   % Uses default DWI group tracts, plots lesion volumes, separate hemispheres
    %   f3('norm', 'hemi')     % Uses default DWI group tracts, plots normalized volumes, separate hemispheres
    %
    % This script creates polar plots:
    % - With 'hemi' option: Two polar plots (left and right hemispheres)
    % - Without 'hemi' option: One polar plot (combined hemispheres)
    %
    % Groups plotted:
    % - Association
    % - Cerebellar
    % - Occipitoparietal
    % - ProjectionBrainstem (PB)
    
    % Parse input arguments
    p = inputParser;
    
    % Handle different input combinations
    if nargin == 1 && ismember(varargin{1}, {'number', 'volume', 'norm'})
        tractType = 'dwi';  % Use default tract type
        metricType = varargin{1};
        separateHemi = false;
    elseif nargin == 2 && ismember(varargin{1}, {'number', 'volume', 'norm'}) && strcmp(varargin{2}, 'hemi')
        tractType = 'dwi';  % Use default tract type
        metricType = varargin{1};
        separateHemi = true;
    else
        % Normal case: parse all parameters
        addOptional(p, 'tractType', 'dwi', @(x) ismember(x, {'dwi', 'template'}));
        addOptional(p, 'metricType', 'number', @(x) ismember(x, {'number', 'volume', 'norm'}));
        addOptional(p, 'hemi', '', @(x) strcmp(x, 'hemi'));
        parse(p, varargin{:});
        tractType = p.Results.tractType;
        metricType = p.Results.metricType;
        separateHemi = strcmp(p.Results.hemi, 'hemi');
    end

    % Add helper functions to path
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));
    
    % Load configuration
    config = load_config();
    
    % Get the correct Excel file path from config
    excelFile = fullfile(config.bidsDir, config.statsDir, ...
                        config.excelOutputs.(tractType).group);
    
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
        else
            warning('Unknown MS type for subject: %s', subID);
            msTypes{i} = 'Unknown';
        end
    end
    data.MSType = categorical(msTypes);
    
    % Define the groups and their display names
    groups = {'Association', 'Cerebellar', 'Occipitoparietal', 'ProjectionBrainstem'};
    groupDisplayNames = {'Association', 'Cerebellar', 'Occipitoparietal', 'ProjectionBrainstem'};
    
    % Initialize arrays for metrics
    numGroups = length(groups);
    if separateHemi
        metricsL = zeros(numGroups, 3); % 3 MS types
        metricsR = zeros(numGroups, 3);
        seL = zeros(numGroups, 3);
        seR = zeros(numGroups, 3);
    else
        metrics = zeros(numGroups, 3);
        se = zeros(numGroups, 3);
    end
    
    % Define metric suffix based on type
    switch metricType
        case 'number'
            suffix = 'LN';
            titleMetric = 'Lesion Numbers';
        case 'volume'
            suffix = 'LV';
            titleMetric = 'Lesion Volumes';
        case 'norm'
            suffix = 'Lnorm';
            titleMetric = 'Normalized Lesion Volumes';
    end
    
    % Calculate averages and standard errors for each group and MS type
    msTypesList = {'RRMS', 'SPMS', 'PPMS'};
    for i = 1:numGroups
        for j = 1:length(msTypesList)
            msType = msTypesList{j};
            idx = data.MSType == msType;
            
            % Get column name
            colName = groups{i};
            if strcmp(colName, 'ProjectionBrainstem')
                colName = 'PB';
            end
            
            if separateHemi
                % Get left hemisphere metrics
                colL = [colName 'L' suffix];
                if ismember(colL, data.Properties.VariableNames)
                    values = data{idx, colL};
                    metricsL(i,j) = mean(values, 'omitnan');
                    seL(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
                
                % Get right hemisphere metrics
                colR = [colName 'R' suffix];
                if ismember(colR, data.Properties.VariableNames)
                    values = data{idx, colR};
                    metricsR(i,j) = mean(values, 'omitnan');
                    seR(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            else
                % Get combined hemisphere metrics
                col = [colName suffix];
                if ismember(col, data.Properties.VariableNames)
                    values = data{idx, col};
                    metrics(i,j) = mean(values, 'omitnan');
                    se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
                end
            end
        end
    end
    
    % Create the figure
    if separateHemi
        fig = figure('Position', [100 100 1200 500]);
        maxCommonR = max(max([metricsL(:) + seL(:); metricsR(:) + seR(:)])) * 1.1;
    else
        fig = figure('Position', [100 100 600 500]);
        maxCommonR = max(max(metrics(:) + se(:))) * 1.1;
    end
    
    if maxCommonR == 0
        maxCommonR = 1; % Handle case of all zeros
    end
    
    if separateHemi
        % Plot left hemisphere
        subplot(1,2,1);
        createPolarPlot(metricsL, seL, groupDisplayNames, msTypesList, config, ...
                        ['Left Hemisphere - ' titleMetric], tractType, maxCommonR, false);
        
        % Plot right hemisphere
        subplot(1,2,2);
        createPolarPlot(metricsR, seR, groupDisplayNames, msTypesList, config, ...
                        ['Right Hemisphere - ' titleMetric], tractType, maxCommonR, false);
        
        % Add a common title
        if strcmp(tractType, 'template')
            titleStr = sprintf('Tract Group %s by MS Type and Hemisphere (Atlas)', titleMetric);
        else
            titleStr = sprintf('Tract Group %s by MS Type and Hemisphere', titleMetric);
        end
    else
        % Plot combined hemispheres
        createPolarPlot(metrics, se, groupDisplayNames, msTypesList, config, ...
                        titleMetric, tractType, maxCommonR, false);
        
        % Add title
        if strcmp(tractType, 'template')
            titleStr = sprintf('Tract Group %s by MS Type (Atlas)', titleMetric);
        else
            titleStr = sprintf('Tract Group %s by MS Type', titleMetric);
        end
    end
    sgtitle(titleStr, 'FontSize', 14);
end
    
  