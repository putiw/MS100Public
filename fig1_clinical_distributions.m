function fig1_clinical_distributions()
clearvars; clc; close all;

   % Add helper functions to path
   addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'helpers')));
    
   % Load configuration
   config = load_config();
   
   % Define clinical scores of interest
   clinical_scores = {'EDSS', 'T25FW', 'x9HPTD', 'x9HPTND', 'SDMTcorrect', 'MSPro'};
   
   % Load clinical data
   clinical_file = fullfile(config.bidsDir, config.statsDir, 'clinicalScore.xlsx');
   clinical_data = readtable(clinical_file);
   
   % Extract MS types from subject IDs
   msTypes = cell(height(clinical_data), 1);
   for i = 1:height(clinical_data)
       subID = clinical_data.SubjectID{i};
       if contains(subID, 'PP')
           msTypes{i} = 'PPMS';
       elseif contains(subID, 'RR')
           msTypes{i} = 'RRMS';
       elseif contains(subID, 'SP')
           msTypes{i} = 'SPMS';
       elseif contains(subID, 'C')
           msTypes{i} = 'control';
       else
           warning('Unknown MS type for subject: %s', subID);
           msTypes{i} = 'Unknown';
       end
   end
   clinical_data.MSType = categorical(msTypes);
   
   % Set default figure properties
   set(0, 'DefaultAxesLineWidth', 2);
   set(0, 'DefaultAxesTickDir', 'out');
   set(0, 'DefaultAxesFontSize', 15);
   set(0, 'DefaultAxesFontName', 'Arial');
   
   % plot the distribution of each the clinical scores
   figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);
   for i = 1:length(clinical_scores)
       subplot(2, 3, i);
       
       % Get data for each MS type in desired order
       msTypesList = {'control', 'RRMS', 'PPMS', 'SPMS'};
       
       % Create histogram bins
       [counts, edges] = histcounts(clinical_data.(clinical_scores{i}));
       
       % Initialize matrix for stacked data
       stacked_data = zeros(length(counts), length(msTypesList));
       
       % Calculate counts for each MS type in each bin
       for j = 1:length(msTypesList)
           msType = msTypesList{j};
           idx = clinical_data.MSType == msType;
           [type_counts, ~] = histcounts(clinical_data.(clinical_scores{i})(idx), edges);
           stacked_data(:,j) = type_counts;
       end
       
       % Create stacked bar plot
       b = bar(edges(1:end-1), stacked_data, 'stacked');
       
       % Set colors for each MS type
       for j = 1:length(msTypesList)
           msType = msTypesList{j};
           b(j).FaceColor = config.msTypeColors.(msType) / 255;
       end
       
       title(clinical_scores{i}, 'FontSize', 15, 'FontWeight', 'bold');
       
       % Set appropriate x-axis label based on the measure
       if ismember(clinical_scores{i}, {'T25FW', 'x9HPTD', 'x9HPTND'})
           xlabel('Time (seconds)', 'FontSize', 15, 'FontAngle', 'italic');
       else
           xlabel('Score', 'FontSize', 15, 'FontAngle', 'italic');
       end
       
       ylabel('Count', 'FontSize', 15, 'FontAngle', 'italic');
       
       % Set axis properties
       ax = gca;
       ax.LineWidth = 2;
       ax.TickDir = 'out';
       ax.FontSize = 15;
       ax.TickLabelInterpreter = 'tex';
       ax.XTickLabel = cellfun(@(x) ['\it' x], ax.XTickLabel, 'UniformOutput', false);
       %ax.YTickLabel = cellfun(@(x) ['\it' x], ax.YTickLabel, 'UniformOutput', false);
       
       % Fix y-axis scaling for each subplot
       max_count = max(sum(stacked_data, 2));
       % Round up to nearest multiple of 10
       max_count = ceil(max_count/10) * 10;
       yticks(0:10:max_count);
       ylim([0 max_count]);
       
       legend(msTypesList, 'Location', 'best', 'FontSize', 12);
   end

% Create radar plot for group differences
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

% Define the order of clinical scores for the radar plot
radar_scores = {'x9HPTND','DMT', 'MSPro', 'x9HPTD', 'EDSS', 'T25FW'};
msTypesList = {'RRMS', 'PPMS', 'SPMS'};

% Calculate mean values and standard errors for each measure
mean_scores = zeros(length(msTypesList), length(radar_scores));
se_scores = zeros(length(msTypesList), length(radar_scores));

% Define hardcoded y-limits for each measure
ylims = [
    0, 40;    
    30, 80;   
    0, 3;   
    0, 40;   
    0, 7;    
    0, 20   
];

for j = 1:length(radar_scores)
    % Get all values for this measure
    if strcmp(radar_scores{j}, 'DMT')
        % Convert SDMT to DMT (incorrect questions)
        all_values = 110 - clinical_data.SDMTcorrect;
    else
        all_values = clinical_data.(radar_scores{j});
    end
    
    % Calculate mean and SE for each group
    for i = 1:length(msTypesList)
        msType = msTypesList{i};
        idx = clinical_data.MSType == msType;
        if strcmp(radar_scores{j}, 'DMT')
            group_values = 110 - clinical_data.SDMTcorrect(idx);
        else
            group_values = clinical_data.(radar_scores{j})(idx);
        end
        mean_scores(i,j) = mean(group_values, 'omitnan');
        se_scores(i,j) = std(group_values, 'omitnan') / sqrt(sum(~isnan(group_values)));
    end
end

% Create polar plot using the helper function
createPolarPlot_diff(mean_scores', se_scores', radar_scores, msTypesList, config, ...
                    'Clinical Scores by MS Type', 'dwi', ylims, false);

