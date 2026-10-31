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
    'Supplementary_FigureS1.pdf');
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
% The original six panels and cell values, at their final publication size.
sizeIn = [169.8 209.8]/25.4;
fig = figure('Visible','off','Color','w','Units','inches', ...
    'Position',[.5 .5 sizeIn],'PaperUnits','inches','PaperSize',sizeIn, ...
    'PaperPosition',[0 0 sizeIn],'DefaultAxesFontName','Arial', ...
    'DefaultTextFontName','Arial');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
palette=flipud(hot(256));
colormap(fig,palette);
set(fig,'Renderer','painters');
for pi = 1:6
    panel = panels(pi);
    col = mod(pi-1,2); row = floor((pi-1)/2);
    ax = axes(fig,'Position',[.085+.50*col .745-.295*row .335 .185]);
    hold(ax,'on');
    [nr,nc]=size(panel.corr_matrix);
    for ii=1:nr
        for jj=1:nc
            v=panel.corr_matrix(ii,jj);
            index=min(256,max(1,floor(v*256)+1));
            patch(ax,jj+[-.5 .5 .5 -.5],ii+[-.5 -.5 .5 .5],palette(index,:), ...
                'EdgeColor','none');
        end
    end
    set(ax,'YDir','reverse','XLim',[.5 nc+.5],'YLim',[.5 nr+.5]);
    clim(ax,[0 1]); box(ax,'on');
    names = panel.display_names;
    for k=1:numel(names)
        names{k}=strrep(names{k},'Association','Assoc.');
        names{k}=strrep(names{k},'Cerebellar','Cereb.');
        names{k}=strrep(names{k},'Occipitoparietal','Occip.-par.');
        names{k}=strrep(names{k},'ProjectionBrainstem','Proj.-BS');
        names{k}=strrep(names{k},'Projection Brainstem','Proj.-BS');
        names{k}=strrep(names{k},'Periventricular','PV');
        names{k}=strrep(names{k},'Juxtacortical','JC');
        names{k}=strrep(names{k},'Infratentorial','IT');
        names{k}=strrep(names{k},'Deep White Matter','DWM');
    end
    set(ax,'XTick',1:numel(panel.clinical_labels),'XTickLabel',panel.clinical_labels, ...
        'YTick',1:numel(names),'YTickLabel',names,'FontSize',7.2, ...
        'TickLabelInterpreter','none','TickLength',[0 0]);
    xtickangle(ax,55);
    title(ax,panel.title,'FontSize',8,'Interpreter','none');
    % Vector colour scale: avoid raster heatmaps and image-based colourbars.
    cb=axes(fig,'Position',[.427+.50*col .745-.295*row .012 .185]);
    hold(cb,'on');
    for kk=1:256
        patch(cb,[0 1 1 0],[(kk-1) (kk-1) kk kk]/256, ...
            palette(kk,:),'EdgeColor','none');
    end
    clim(cb,[0 1]);
    set(cb,'XLim',[0 1],'YLim',[0 1],'XTick',[], ...
        'YTick',[0 .5 1],'YAxisLocation','right','FontName','Arial', ...
        'FontSize',7.2,'TickLength',[0 0],'Box','on');
    ylabel(cb,'|partial r|','FontSize',7.2);
    for i=1:size(panel.corr_matrix,1)
        for j=1:size(panel.corr_matrix,2)
            v=panel.corr_matrix(i,j); if isnan(v),continue;end
            label=sprintf('%.2f',v); weight='normal';
            if panel.significant(i,j),label=[label '*'];weight='bold';end
            color='k';if v>.6,color='w';end
            text(ax,j,i,label,'HorizontalAlignment','center','FontSize',7.2, ...
                'FontWeight',weight,'Color',color);
        end
    end
end
annotation(fig,'textbox',[.035 .01 .93 .055],'String', ...
    {'Assoc.: Association; Cereb.: Cerebellar; Occip.-par.: Occipitoparietal; Proj.-BS: Projection/Brainstem', ...
    'PV: Periventricular; JC: Juxtacortical; IT: Infratentorial; DWM: Deep white matter', ...
    '* Bonferroni significant (182 tests). Colours and numbers: absolute partial correlation (unitless).'}, ...
    'FontName','Arial','FontSize',7.2,'EdgeColor','none','Interpreter','none');
print(fig,pdfPath,'-dpdf','-painters');
print(fig,pngPath,'-dpng','-r300');
end
