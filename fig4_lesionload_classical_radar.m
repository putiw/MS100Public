function fig4_lesionload_classical_radar(varargin)
    % F3_ICO Creates polar plots of icometrix lesion metrics by MS type
    %
    % Usage:
    %   f3_ico()                    % Uses default lesion numbers
    %   f3_ico('number')            % Plots lesion numbers
    %   f3_ico('volume')            % Plots lesion volumes
    %   f3_ico('norm')
    %
    % This script creates a polar plot showing:
    % Average metric per region for each MS type
    %
    % Regions plotted:
    % - Periventricular
    % - Juxtacortical
    % - Infratentorial
    % - Deep White Matter
    
    % Parse input arguments
    p = inputParser;
    addOptional(p, 'metricType', 'number', @(x) ismember(x, {'number', 'volume', 'norm'}));
    parse(p, varargin{:});
    metricType = p.Results.metricType;

    % Add helper functions to path
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));
    
    % Load configuration
    config = load_config();
    
    % Get the Excel file path
    excelFile = fullfile(config.bidsDir, config.statsDir, 'icometrixLesionLoad.xlsx');
    
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
    
    % Define the regions
    regions = {'periventricular', 'juxtacortical', 'infratentorial', 'deepwhitematter'};
    regionDisplayNames = {'Periventricular', 'Juxtacortical', 'Infratentorial', 'Deep WM'};
    
    % Initialize arrays for metrics
    numRegions = length(regions);
    metrics = zeros(numRegions, 3); % 3 MS types
    
    % Initialize arrays for standard errors
    se = zeros(numRegions, 3);
    
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
    
    % Calculate averages and standard errors for each region and MS type
    msTypesList = {'RRMS', 'SPMS', 'PPMS'};
    for i = 1:numRegions
        for j = 1:length(msTypesList)
            msType = msTypesList{j};
            idx = data.MSType == msType;
            
            % Get metrics
            colName = [regions{i} suffix];
            if ismember(colName, data.Properties.VariableNames)
                values = data{idx, colName};
                metrics(i,j) = mean(values, 'omitnan');
                se(i,j) = std(values, 'omitnan') / sqrt(sum(~isnan(values)));
            end
        end
    end
    
    % Create the figure
    fig = figure('Position', [100 100 600 500]);
    
    % Calculate max radius for the plot
    maxR = max(max(metrics(:) + se(:))) * 1.1;
    if maxR == 0
        maxR = 1; % Handle case of all zeros
    end
    
    % Create polar plot
    createPolarPlot(metrics, se, regionDisplayNames, msTypesList, config, ...
                    titleMetric, 'icometrix', maxR, false);
    
    % Add title
    title(sprintf('Region %s by MS Type', titleMetric), 'FontSize', 14);
end 