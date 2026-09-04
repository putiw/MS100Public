function report = figure6_r2_comparison_sdmt_msfc_noddi(varargin)
% FIGURE6_R2_COMPARISON_SDMT_MSFC_NODDI Extend manuscript Figure 6/S2.
%
% Copied from MS100Public/rev_fig6_r2_comparison.m. Its grouped bars,
% colors, sample-size labels, matched-subject markers, and typography are
% retained; the explicit analysis-facing changes are limited to:
%   1. five rather than three continuous clinical outcomes;
%   2. seven rather than four imaging metrics;
%   3. a second pass for the age/gender-adjusted supplementary figure; and
%   4. configurable input/output paths for the reviewer package.
% Figure 6 and Supplementary Figure S2 use the final one-row-by-five-column
% layout at the journal's 170-mm page width.
%
% NDI, ODI, and ISOVF use the same tract-based/classical-region comparison as
% the original four metrics, including separate exact matched-subject refits.

p = inputParser;
scriptDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(scriptDir);
addParameter(p, 'fullFile', fullfile(packageRoot, 'results', 'model_results', ...
    'figure6_full_models.csv'), @(x) ischar(x) || isstring(x));
addParameter(p, 'sameNFile', fullfile(packageRoot, 'results', 'model_results', ...
    'figure6_matched_models.csv'), @(x) ischar(x) || isstring(x));
addParameter(p, 'outputDir', fullfile(packageRoot, 'results', 'figures'), ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'previewDir', fullfile(packageRoot, 'work', ...
    'figure_previews'), @(x) ischar(x) || isstring(x));
addParameter(p, 'doExport', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

fullFile = char(p.Results.fullFile);
sameNFile = char(p.Results.sameNFile);
outputDir = char(p.Results.outputDir);
previewDir = char(p.Results.previewDir);
assert(isfile(fullFile), 'Missing: %s', fullFile);
assert(isfile(sameNFile), 'Missing: %s', sameNFile);

% Original lists extended by exactly three metrics and two outcomes.
metricsToShow = {'T1','MTR','FA','MD','NDI','ODI','ISOVF'};
clinicalLabels = {'T25FW','9HPT-D','9HPT-ND','SDMT','MSFC-SDMT'};
clinicalKeys = {'T25FW','x9HPTD','x9HPTND','SDMTcorrect','MSFC_SDMT'};

colorTract = [73 167 246] / 255;
colorClassic = [255 199 119] / 255;
Tfull = readtable(fullFile, 'TextType','string');
TsameN = readtable(sameNFile, 'TextType','string');
assert(height(Tfull) == 70, ...
    'Expected 70 full-model rows (5 outcomes x 7 metrics x 2 models).');
assert(height(TsameN) == 70, ...
    'Expected 70 matched rows (5 outcomes x 7 metrics x 2 models).');
Dfull = parse_table(Tfull);
DsameN = parse_table(TsameN);

if p.Results.doExport && ~isfolder(outputDir), mkdir(outputDir); end
if p.Results.doExport && ~isfolder(previewDir), mkdir(previewDir); end
outputs = render_figures(Dfull,DsameN,metricsToShow,clinicalLabels, ...
    clinicalKeys,colorTract,colorClassic,outputDir,previewDir, ...
    p.Results.doExport);

report = struct();
report.source_script = 'code/figure6_r2_comparison_sdmt_msfc_noddi.m';
report.metrics = metricsToShow;
report.outcomes = clinicalLabels;
report.figure6 = outputs.figure6;
report.supplementary_figure_s2 = outputs.supplementary_figure_s2;
report.layout = outputs.layout;
report.noddi_framework = 'tract-based and classical-region';
report.exported = p.Results.doExport;
end


function D = parse_table(T)
% Original Figure 6 long-table parser, with R2_Demo retained for S2.
vn = T.Properties.VariableNames;
col = @(c) vn{strcmpi(vn,c)};
D.key = string(T.(col('ClinicalMetric')));
D.model = lower(strtrim(string(T.(col('Model')))));
D.N = to_num(T.(col('N')));
D.R2 = to_num(T.(col('R2')));
D.R2_Demo = to_num(T.(col('R2_Demo')));
end


function v = to_num(v)
if isnumeric(v), v = double(v); return; end
if iscell(v), v = string(v); end
if isstring(v), v = str2double(v); return; end
v = double(v);
end


function outputs = render_figures(Dfull,DsameN,metrics,clinicalLabels, ...
        clinicalKeys,colorTract,colorClassic,outputDir,previewDir,doExport)
% Final vertical grouped-bar layout, using the same frozen model tables as
% the original manuscript figure.

valueColumns = {'R2','R2_Demo'};
titles = {'R^2 Comparison: Tract-based vs Classic (without demographics)', ...
          'R^2 Comparison: Tract-based vs Classic (with age and recorded gender)'};
pdfNames = {'Figure6_with_SDMT_MSFC_NODDI.pdf', ...
    'Supplementary_Figure_S2_with_SDMT_MSFC_NODDI.pdf'};
pngNames = {'Figure6_with_SDMT_MSFC_NODDI.png', ...
    'Supplementary_Figure_S2_with_SDMT_MSFC_NODDI.png'};
% 169 x 80 mm: full journal width, inside the 170 x 210 mm limit.
figureSize = [6.65 3.15];
barWidth = 0.92;

for fi = 1:numel(valueColumns)
    valueColumn = valueColumns{fi};
    fig = figure('Visible','off', 'Units','inches', ...
        'Position',[0.5 0.5 figureSize], 'Color','white');
    tl = tiledlayout(fig,1,5,'TileSpacing','compact','Padding','compact');
    title(tl,titles{fi},'FontName','Arial','FontSize',9.5, ...
        'FontWeight','bold','Interpreter','tex');
    firstAx = [];
    legendHandles = gobjects(1,3);

    for ci = 1:numel(clinicalKeys)
        values = collect_panel_values(Dfull,DsameN,clinicalKeys{ci}, ...
            metrics,valueColumn);
        ax = nexttile(tl,ci);
        hold(ax,'on');
        b = bar(ax,[values.r2Tract(:),values.r2Classic(:)], ...
            'grouped','BarWidth',barWidth);
        b(1).FaceColor = colorTract;
        b(2).FaceColor = colorClassic;
        b(1).EdgeColor = 'none';
        b(2).EdgeColor = 'none';
        drawnow;

        for mi = 1:numel(metrics)
            xT = b(1).XEndPoints(mi);
            xC = b(2).XEndPoints(mi);
            topT = max([values.r2Tract(mi),values.r2TractSameN(mi)], ...
                [],'omitnan');
            topC = max([values.r2Classic(mi),values.r2ClassicSameN(mi)], ...
                [],'omitnan');
            add_one_by_five_label(ax,xT,topT,values.nTract(mi));
            add_one_by_five_label(ax,xC,topC,values.nClassic(mi));
            add_one_by_five_marker(ax,xT,values.r2TractSameN(mi));
            add_one_by_five_marker(ax,xC,values.r2ClassicSameN(mi));
        end

        hold(ax,'off');
        box(ax,'off');
        grid(ax,'off');
        set(ax,'FontName','Arial','LineWidth',0.85,'FontSize',7, ...
            'FontAngle','italic','TickDir','out');
        title(ax,clinicalLabels{ci},'FontName','Arial','FontSize',8.5, ...
            'FontWeight','bold','FontAngle','italic','Interpreter','none');
        xticks(ax,1:numel(metrics));
        xticklabels(ax,metrics);
        xtickangle(ax,90);
        xlim(ax,[0.45,numel(metrics)+0.55]);
        ylim(ax,[0,0.56]);
        yticks(ax,0:0.1:0.5);
        if ci == 1
            ylabel(ax,'R^2', ...
                'FontName','Arial','FontSize',8.5, ...
                'FontWeight','bold','FontAngle','italic','Interpreter','tex');
            firstAx = ax;
            legendHandles(1:2) = b;
            hold(ax,'on');
            legendHandles(3) = plot(ax,nan,nan,'k-','LineWidth',2.2);
            hold(ax,'off');
        else
            yticklabels(ax,{});
        end
    end

    lg = legend(firstAx,legendHandles, ...
        {'Tract-based','Classic','Matched N'}, ...
        'Orientation','horizontal','Location','northoutside', ...
        'FontName','Arial','FontSize',7.5,'FontAngle','italic','Box','off');
    lg.Layout.Tile = 'north';
    drawnow;
    set(fig,'Position',[0.5 0.5 figureSize]);
    drawnow;
    if doExport
        set(fig,'PaperUnits','inches','PaperSize',figureSize, ...
            'PaperPosition',[0 0 figureSize],'PaperPositionMode','manual');
        set(fig,'InvertHardcopy','off');
        print(fig,fullfile(outputDir,pdfNames{fi}),'-dpdf','-painters');
        drawnow;
        exportgraphics(fig,fullfile(previewDir,pngNames{fi}), ...
            'Resolution',300,'BackgroundColor','white');
    end
    close(fig);
end

outputs = struct( ...
    'figure6',fullfile(outputDir,pdfNames{1}), ...
    'supplementary_figure_s2',fullfile(outputDir,pdfNames{2}), ...
    'layout','one row by five columns', ...
    'bar_width',barWidth, ...
    'exported',doExport);
end


function values = collect_panel_values(Dfull,DsameN,clinicalKey,metrics,valueColumn)
nMetrics = numel(metrics);
values.r2Tract = nan(1,nMetrics);
values.r2Classic = nan(1,nMetrics);
values.nTract = nan(1,nMetrics);
values.nClassic = nan(1,nMetrics);
values.r2TractSameN = nan(1,nMetrics);
values.r2ClassicSameN = nan(1,nMetrics);
for mi = 1:nMetrics
    key = sprintf('%s-%s',clinicalKey,metrics{mi});
    rowT = find(strcmp(Dfull.key,key) & strcmp(Dfull.model,'tract-based'));
    rowC = find(strcmp(Dfull.key,key) & strcmp(Dfull.model,'classic'));
    rowTs = find(strcmp(DsameN.key,key) & strcmp(DsameN.model,'tract-based'));
    rowCs = find(strcmp(DsameN.key,key) & strcmp(DsameN.model,'classic'));
    assert(isscalar(rowT) && isscalar(rowC), ...
        'Expected one full tract/classical pair for %s.',key);
    assert(isscalar(rowTs) && isscalar(rowCs), ...
        'Expected one matched tract/classical pair for %s.',key);
    assert(DsameN.N(rowTs) == DsameN.N(rowCs), ...
        'Matched tract/classical sample sizes differ for %s.',key);
    values.r2Tract(mi) = Dfull.(valueColumn)(rowT);
    values.r2Classic(mi) = Dfull.(valueColumn)(rowC);
    values.nTract(mi) = Dfull.N(rowT);
    values.nClassic(mi) = Dfull.N(rowC);
    values.r2TractSameN(mi) = DsameN.(valueColumn)(rowTs);
    values.r2ClassicSameN(mi) = DsameN.(valueColumn)(rowCs);
end
end


function add_one_by_five_label(ax,xValue,topValue,nValue)
if isfinite(topValue) && isfinite(nValue)
    text(ax,xValue,topValue+0.010,sprintf('%d',round(nValue)), ...
        'HorizontalAlignment','center','FontName','Arial','FontSize',6.5, ...
        'FontAngle','italic','Color','k');
end
end


function add_one_by_five_marker(ax,xValue,r2Value)
if isfinite(r2Value)
    plot(ax,[xValue-0.11,xValue+0.11],[r2Value,r2Value], ...
        'k-','LineWidth',1.3);
end
end
