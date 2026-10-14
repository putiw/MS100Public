function [barFig, polarFig, summary] = fig4_plain_matched89(outputDir, tractFile, classicalFile)
%FIG4_PLAIN_MATCHED89 Reproduce the numerical Figure 4 as two MATLAB figures.
%
% In the original manuscript, the classical lesion-load workbook contained
% 50 RRMS participants, but the tract workbook contained only 49 RRMS
% participants. The extra classical participant (sub-RR047) had no usable
% diffusion-derived tract data. The original classical and tract plots
% therefore did not use exactly the same people, despite the Figure 4
% caption reporting 49 RRMS. This revised Figure 4 uses the 89 participants
% with tract data for BOTH frameworks: RRMS 49, PPMS 17, SPMS 23.
%
% Output 1 (Figure 4A): grouped whole-brain and four-tract-region lesion
% volume bars, mean +/- sample SD, in mL, with a 0-50 mL y-axis.
% Output 2 (Figure 4B): six polar panels (classical/tract x LN/LV/Lnorm),
% regional mean +/- SEM. Zeros are retained. Each radial maximum is
% 1.1 * max(mean + SEM), as in the original MATLAB radar scripts.
%
% These are plain numerical plots, not the final Illustrator composition.
% The four tract-region volumes are summed to reproduce the original bar;
% overlapping tract masks may contribute to more than one region.
%
% Usage:
%   fig4_plain_matched89()                      % Open two MATLAB figures
%   fig4_plain_matched89('/path/to/new/output') % Save two vector PDFs
%   fig4_plain_matched89(out, tractXlsx, classicalXlsx)
%
% The source spreadsheets are read only. The export directory must not
% exist, so an earlier figure can never be silently overwritten.

if nargin < 1
    outputDir = '';
end
repoDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
config = jsondecode(fileread(fullfile(repoDir, 'configs', 'config.json')));
statsDir = fullfile(config.bidsDir, config.statsDir);
if nargin < 2 || isempty(tractFile)
    tractFile = fullfile(statsDir, config.excelOutputs.dwi.group);
end
if nargin < 3 || isempty(classicalFile)
    classicalFile = fullfile(statsDir, 'icometrixLesionLoad.xlsx');
end
if ~isfile(tractFile) || ~isfile(classicalFile)
    error('Both Figure 4 source workbooks are required: %s ; %s', tractFile, classicalFile);
end
if ~isempty(outputDir) && (isfolder(outputDir) || isfile(outputDir))
    error('Refusing to overwrite an existing output path: %s', outputDir);
end

tract = readtable(tractFile);
classical = readtable(classicalFile);
if ~ismember('SubjectID', tract.Properties.VariableNames) || ...
        ~ismember('SubjectID', classical.Properties.VariableNames)
    error('Both workbooks must have a SubjectID column');
end
tractIDs = string(tract.SubjectID);
classicalIDs = string(classical.SubjectID);
if any(ismissing(tractIDs)) || any(ismissing(classicalIDs)) || ...
        numel(unique(tractIDs)) ~= 89 || numel(unique(classicalIDs)) ~= 90
    error('Expected 89 unique tract IDs and 90 unique classical IDs');
end
[found, classicalRows] = ismember(tractIDs, classicalIDs);
if ~all(found) || ~isequal(setdiff(classicalIDs, tractIDs), "sub-RR047")
    error('Figure 4 cohort mismatch: RR047 must be the sole classical-only participant');
end
classical = classical(classicalRows, :);

groups = {'RRMS', 'PPMS', 'SPMS'};
prefixes = ["sub-RR", "sub-PP", "sub-SP"];
expectedN = [49, 17, 23];
n = zeros(1, 3);
groupMask = false(height(tract), 3);
for g = 1:3
    groupMask(:, g) = startsWith(tractIDs, prefixes(g));
    n(g) = sum(groupMask(:, g));
end
if ~isequal(n, expectedN) || any(sum(groupMask, 2) ~= 1)
    error('Matched Figure 4 counts are not RRMS 49, PPMS 17, SPMS 23');
end

classicalRegions = {'periventricular', 'juxtacortical', 'infratentorial', 'deepwhitematter'};
tractRegions = {'Association', 'Cerebellar', 'Occipitoparietal', 'PB'};
regions = {classicalRegions, tractRegions};
regionLabels = {{'Periventricular', 'Juxtacortical', 'Infratentorial', 'Deep WM'}, ...
                {'Association', 'Cerebellar', 'Occipitoparietal', 'Projection/Brainstem'}};
metricSuffix = {'LN', 'LV', 'Lnorm'};
metricLabels = {'Lesion number', 'Lesion volume (mL)', 'Normalized lesion volume'};
sources = {classical, tract};

wholeBrain = zeros(height(tract), 1);
tractRegionSum = zeros(height(tract), 1);
for r = 1:4
    wholeBrain = wholeBrain + checkedColumn(classical, [classicalRegions{r}, 'LV']);
    tractRegionSum = tractRegionSum + checkedColumn(tract, [tractRegions{r}, 'LV']);
end
if any(abs(wholeBrain - checkedColumn(tract, 'WBLV')) > 1e-8)
    error('The tract WBLV column does not equal the classical whole-brain LV sum');
end
wholeBrain = wholeBrain / 1000;      % mm^3 to mL
tractRegionSum = tractRegionSum / 1000;
barMeans = zeros(3, 2);
barSD = zeros(3, 2);
for g = 1:3
    barMeans(g, :) = [mean(wholeBrain(groupMask(:, g))), ...
                      mean(tractRegionSum(groupMask(:, g)))];
    barSD(g, :) = [std(wholeBrain(groupMask(:, g)), 0), ...
                   std(tractRegionSum(groupMask(:, g)), 0)];
end
if any(barMeans(:) + barSD(:) > 50)
    error('Figure 4A values exceed the original 0-50 mL plotting range');
end

% Dimensions: region x subtype x framework x metric.
radarMeans = zeros(4, 3, 2, 3);
radarSEM = zeros(4, 3, 2, 3);
radialMax = zeros(2, 3);
for f = 1:2
    for m = 1:3
        for r = 1:4
            column = [regions{f}{r}, metricSuffix{m}];
            values = checkedColumn(sources{f}, column);
            for g = 1:3
                selected = values(groupMask(:, g));
                divisor = 1;
                if strcmp(metricSuffix{m}, 'LV')
                    divisor = 1000;
                end
                selected = selected / divisor;
                radarMeans(r, g, f, m) = mean(selected);
                radarSEM(r, g, f, m) = std(selected, 0) / sqrt(numel(selected));
            end
        end
        cellUpper = radarMeans(:, :, f, m) + radarSEM(:, :, f, m);
        radialMax(f, m) = 1.1 * max(cellUpper(:));
    end
end

summary = struct('groups', {groups}, 'n', n, 'barMeans_mL', barMeans, ...
                 'barSD_mL', barSD, 'radarMeans', radarMeans, ...
                 'radarSEM', radarSEM, 'radialMax', radialMax, ...
                 'classicalRegions', {classicalRegions}, 'tractRegions', {tractRegions}, ...
                 'metrics', {metricSuffix});
palette = [reshape(config.msTypeColors.RRMS, 1, 3); ...
           reshape(config.msTypeColors.PPMS, 1, 3); ...
           reshape(config.msTypeColors.SPMS, 1, 3)] / 255;

visibility = 'on';
if ~isempty(outputDir)
    visibility = 'off';
end
barFig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', ...
                'Position', [1, 1, 8, 5]);
axBar = axes(barFig);
hold(axBar, 'on');
bars = bar(axBar, barMeans, 'grouped');
bars(1).FaceColor = [1, 1, 1];
bars(1).EdgeColor = [0.2, 0.2, 0.2];
bars(2).FaceColor = [0.97, 0.73, 0.41];
bars(2).EdgeColor = [0.65, 0.36, 0.04];
for k = 1:2
    errorbar(axBar, bars(k).XEndPoints, barMeans(:, k), barSD(:, k), ...
             'k', 'LineStyle', 'none', 'CapSize', 8, 'LineWidth', 1);
end
xticks(axBar, 1:3);
groupTickLabels = arrayfun(@(g) sprintf('%s (n=%d)', groups{g}, n(g)), ...
                           1:3, 'UniformOutput', false);
xticklabels(axBar, groupTickLabels);
ylabel(axBar, 'Lesion volume (mL)');
ylim(axBar, [0, 50]);
yticks(axBar, 0:10:50);
title(axBar, 'Figure 4A. Lesion volume (mean \pm SD)');
legend(axBar, bars, {'Whole brain', 'Sum of four tract regions'}, ...
       'Location', 'northwest', 'Box', 'off');
grid(axBar, 'on');
axBar.XGrid = 'off';
axBar.Layer = 'bottom';
box(axBar, 'off');

polarFig = figure('Visible', visibility, 'Color', 'w', 'Units', 'inches', ...
                  'Position', [1, 1, 15, 10]);
layout = tiledlayout(polarFig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, 'Figure 4B. Regional lesion load (mean \pm SEM)');
frameworkLabels = {'Classical', 'Tract'};
for f = 1:2
    for m = 1:3
        ax = nexttile(layout);
        lineHandles = drawRadar(ax, radarMeans(:, :, f, m), ...
                                radarSEM(:, :, f, m), radialMax(f, m), ...
                                regionLabels{f}, palette);
        title(ax, sprintf('%s: %s', frameworkLabels{f}, metricLabels{m}));
        if f == 1 && m == 1
            firstAx = ax;
            firstLines = lineHandles;
        end
    end
end
lgd = legend(firstAx, firstLines, groups, 'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'south';

if ~isempty(outputDir)
    [created, message] = mkdir(outputDir);
    if ~created
        error('Could not create Figure 4 output directory: %s', message);
    end
    exportgraphics(barFig, fullfile(outputDir, 'Figure4A_plain_matched89.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(polarFig, fullfile(outputDir, 'Figure4B_plain_matched89.pdf'), ...
                   'ContentType', 'vector');
end
end

function values = checkedColumn(data, column)
if ~ismember(column, data.Properties.VariableNames)
    error('Required Figure 4 column is missing: %s', column);
end
values = data{:, column};
if ~isnumeric(values) || any(~isfinite(values))
    error('Missing or nonnumeric Figure 4 values in: %s', column);
end
values = double(values(:));
end

function lineHandles = drawRadar(ax, means, sems, maxRadius, labels, palette)
theta = (0:3)' * (pi / 2) + pi / 4;
circleTheta = linspace(0, 2*pi, 100);
hold(ax, 'on');
for ring = 1:5
    radius = ring / 5;
    plot(ax, radius * cos(circleTheta), radius * sin(circleTheta), ...
         '--', 'Color', [0.72, 0.72, 0.72], 'LineWidth', 0.7, ...
         'HandleVisibility', 'off');
    text(ax, radius + 0.02, 0, sprintf('%.1f', maxRadius * radius), ...
         'FontSize', 7, 'HorizontalAlignment', 'left');
end
for r = 1:4
    plot(ax, [0, cos(theta(r))], [0, sin(theta(r))], '--', ...
         'Color', [0.72, 0.72, 0.72], 'LineWidth', 0.7, ...
         'HandleVisibility', 'off');
    x = 1.16 * cos(theta(r));
    y = 1.16 * sin(theta(r));
    horizontal = 'left';
    if x < 0
        horizontal = 'right';
    end
    text(ax, x, y, labels{r}, 'FontSize', 8, ...
         'HorizontalAlignment', horizontal, 'VerticalAlignment', 'middle');
end
lineHandles = gobjects(1, 3);
for g = [1, 3, 2] % Original radar layering: RRMS, SPMS, PPMS.
    center = means(:, g) / maxRadius;
    upper = (means(:, g) + sems(:, g)) / maxRadius;
    lower = max(0, means(:, g) - sems(:, g)) / maxRadius;
    xPatch = [upper .* cos(theta); flipud(lower .* cos(theta))];
    yPatch = [upper .* sin(theta); flipud(lower .* sin(theta))];
    patch(ax, xPatch, yPatch, palette(g, :), 'FaceAlpha', 0.3, ...
          'EdgeColor', 'none', 'HandleVisibility', 'off');
    closed = [center; center(1)];
    closedTheta = [theta; theta(1)];
    lineHandles(g) = plot(ax, closed .* cos(closedTheta), ...
                          closed .* sin(closedTheta), 'Color', palette(g, :), ...
                          'LineWidth', 2);
end
axis(ax, 'equal');
axis(ax, 'off');
xlim(ax, [-1.35, 1.35]);
ylim(ax, [-1.35, 1.35]);
end
