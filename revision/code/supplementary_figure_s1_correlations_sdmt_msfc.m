function report = supplementary_figure_s1_correlations_sdmt_msfc(varargin)
% SUPPLEMENTARY_FIGURE_S1_CORRELATIONS_SDMT_MSFC Reproduce revised S1.
%
% The plotting section is copied from the validated SDMT/MSFC correlation
% workflow. It reads frozen aggregate correlation results; no participant-
% level data are included in this package.

p = inputParser;
scriptDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(scriptDir);
addParameter(p, 'resultsFile', fullfile(packageRoot, 'results', ...
    'model_results', 'supplementary_figure_s1_correlations.json'), ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'outputDir', fullfile(packageRoot, 'results', 'figures'), ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'previewDir', fullfile(packageRoot, 'work', ...
    'figure_previews'), @(x) ischar(x) || isstring(x));
addParameter(p, 'doExport', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

resultsFile = char(p.Results.resultsFile);
outputDir = char(p.Results.outputDir);
previewDir = char(p.Results.previewDir);
assert(isfile(resultsFile), 'Missing: %s', resultsFile);

payload = jsondecode(fileread(resultsFile));
cells = payload.cells;
assert(numel(cells) == 182, 'Expected 182 correlation cells.');
assert(payload.multiple_testing.global_family.multiplier == 182, ...
    'Expected a single 182-test Bonferroni family.');

panels = build_panels(cells);
pdfPath = fullfile(outputDir, ...
    'Supplementary_Figure_S1_with_SDMT_MSFC.pdf');
pngPath = fullfile(previewDir, ...
    'Supplementary_Figure_S1_with_SDMT_MSFC.png');

if p.Results.doExport
    if ~isfolder(outputDir), mkdir(outputDir); end
    if ~isfolder(previewDir), mkdir(previewDir); end
    render_composite(panels, pdfPath, pngPath);
end

report = struct();
report.source_script = ...
    'code/supplementary_figure_s1_correlations_sdmt_msfc.m';
report.results_file = resultsFile;
report.figure = pdfPath;
report.panels = numel(panels);
report.correlation_cells = numel(cells);
report.bonferroni_multiplier = ...
    payload.multiple_testing.global_family.multiplier;
report.exported = p.Results.doExport;
end


function panels = build_panels(cells)
displayed = cells([cells.DisplayedInS1]);
panelOrders = unique([displayed.PanelOrder], 'stable');
assert(numel(panelOrders) == 6, 'Expected six displayed panels.');
panels = repmat(struct('title', '', 'clinical_labels', {{}}, ...
    'display_names', {{}}, 'corr_matrix', [], 'significant', []), ...
    1, numel(panelOrders));

for panelIndex = 1:numel(panelOrders)
    rows = displayed([displayed.PanelOrder] == panelOrders(panelIndex));
    imagingLabels = unique({rows.ImagingLabel}, 'stable');
    outcomeVariables = unique({rows.OutcomeVariable}, 'stable');
    outcomeLabels = cell(1, numel(outcomeVariables));
    corrMatrix = nan(numel(imagingLabels), numel(outcomeVariables));
    significant = false(size(corrMatrix));

    for outcomeIndex = 1:numel(outcomeVariables)
        first = find(strcmp({rows.OutcomeVariable}, ...
            outcomeVariables{outcomeIndex}), 1, 'first');
        outcomeLabels{outcomeIndex} = rows(first).OutcomeLabel;
    end
    for rowIndex = 1:numel(rows)
        imagingIndex = find(strcmp(imagingLabels, ...
            rows(rowIndex).ImagingLabel), 1);
        outcomeIndex = find(strcmp(outcomeVariables, ...
            rows(rowIndex).OutcomeVariable), 1);
        corrMatrix(imagingIndex, outcomeIndex) = ...
            rows(rowIndex).AbsolutePartialR;
        significant(imagingIndex, outcomeIndex) = ...
            rows(rowIndex).Significant_Bonferroni;
    end

    panels(panelIndex).title = rows(1).PanelTitle;
    panels(panelIndex).clinical_labels = outcomeLabels;
    panels(panelIndex).display_names = imagingLabels;
    panels(panelIndex).corr_matrix = corrMatrix;
    panels(panelIndex).significant = significant;
end
end


function render_composite(panels, pdfPath, pngPath)
% Retain the established S1 dimensions, heatmap, typography and spacing.
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 13 11.5], 'PaperUnits', 'inches', ...
    'PaperSize', [13 11.5], 'PaperPosition', [0 0 13 11.5], ...
    'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial');
closeFigure = onCleanup(@() close(fig)); %#ok<NASGU>
layout = tiledlayout(fig, 3, 2, ...
    'TileSpacing', 'compact', 'Padding', 'loose');
colormap(fig, flipud(hot(256)));

for panelIndex = 1:numel(panels)
    panel = panels(panelIndex);
    ax = nexttile(layout, panelIndex);
    imagesc(ax, panel.corr_matrix);
    clim(ax, [0 1]);
    axis(ax, 'tight');
    box(ax, 'on');
    cb = colorbar(ax);
    cb.FontSize = 7;
    cb.Label.String = '|partial r|';
    cb.Label.FontSize = 7;
    title(ax, [panel.title ' Correlations'], ...
        'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
    set(ax, 'XTick', 1:numel(panel.clinical_labels), ...
        'XTickLabel', panel.clinical_labels, ...
        'YTick', 1:numel(panel.display_names), ...
        'YTickLabel', panel.display_names, ...
        'FontSize', 8, 'TickLabelInterpreter', 'none');
    xtickangle(ax, 27);

    for imagingIndex = 1:size(panel.corr_matrix, 1)
        for outcomeIndex = 1:size(panel.corr_matrix, 2)
            r = panel.corr_matrix(imagingIndex, outcomeIndex);
            if isnan(r), continue; end
            cellText = sprintf('%.2f', r);
            if panel.significant(imagingIndex, outcomeIndex)
                cellText = [cellText '*']; %#ok<AGROW>
                fontWeight = 'bold';
            else
                fontWeight = 'normal';
            end
            if r > 0.6, fontColor = 'w'; else, fontColor = 'k'; end
            text(ax, outcomeIndex, imagingIndex, cellText, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 8, 'FontWeight', fontWeight, ...
                'Color', fontColor);
        end
    end
end

exportgraphics(fig, pdfPath, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
exportgraphics(fig, pngPath, 'Resolution', 300, ...
    'BackgroundColor', 'white');
end
