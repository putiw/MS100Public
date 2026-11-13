function radialLimits = render_supplementary_figure_s3(plottedCsv, outputPdf, limitsJson)
% Plain radar plots from the aggregate NODDI values.
if nargin < 3, limitsJson = ''; end
if isfile(outputPdf), error('Refusing to overwrite %s', outputPdf); end
if ~isempty(limitsJson) && isfile(limitsJson)
    error('Refusing to overwrite %s', limitsJson);
end
D = readtable(plottedCsv, 'TextType', 'string');
metrics = ["NDI", "ODI", "ISOVF"];
allTissues = ["Lesion", "Tail", "NAWM"];
tissues = allTissues(ismember(allTissues, D.Tissue));
groups = ["RRMS", "SPMS", "PPMS", "Control"];
frameworks = ["Classical", "Tract", "Tract"];
tissues = ["Lesion", "Lesion", "NAWM"];
rowKeys = ["ClassicalLesion", "TractLesion", "TractNAWM"];
tractRegions = ["Association", "Cerebellar", "Occipitoparietal", "ProjectionBrainstem"];
regions = {["periventricular", "juxtacortical", "infratentorial", "deepwhitematter"], ...
    tractRegions, tractRegions};
tractLabels = {'Association', 'Cerebellar', {'Occipito','parietal'}, {'Brainstem','Projection'}};
labels = {{'Periventricular', 'Juxtacortical', 'Infratentorial', 'Deep WM'}, ...
    tractLabels, tractLabels};
radialLimits = struct();
for row = 1:numel(tissues)
    for mi = 1:numel(metrics)
        selected = D.Framework == frameworks(row) & D.Tissue == tissues(row) & D.Metric == metrics(mi);
        assert(any(selected),'Missing Figure S3 panel.');
        low = min(D.Mean(selected) - D.SEM(selected));
        high = max(D.Mean(selected) + D.SEM(selected));
        assert(isfinite(low) && isfinite(high) && low < high,'Invalid radial range.');
        radialLimits.(rowKeys(row)).(metrics(mi)) = [floor(100*low), ceil(100*high)]/100;
    end
end

% Use one page coordinate system so the legend cannot resize a radar.
width = 170*72/25.4;
height = 81 + 110*numel(tissues);
sizeIn = [width height]/72;
fig = figure('Visible','off','Color','w','Units','inches','Position',[.5 .5 sizeIn], ...
    'PaperUnits','inches','PaperSize',sizeIn,'PaperPosition',[0 0 sizeIn]);
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
ax = axes(fig,'Position',[0 0 1 1],'XLim',[0 width],'YLim',[0 height], ...
    'DataAspectRatio',[1 1 1],'Visible','off');
hold(ax,'on');
colors = ax.ColorOrder;
theta = (0:4)*pi/2+pi/4;
circle = linspace(0,2*pi,180);
radius = 45;
centersX = [104 253 402];
text(ax,10,height-16,'Classical versus Tract-based NODDI Metrics across MS Subtypes', ...
    'Interpreter','none','VerticalAlignment','top');

for row = 1:numel(tissues)
    tissue = tissues(row);
    framework = frameworks(row);
    cy = height-87-(row-1)*110;
    rowGroups = groups;
    if tissue == "Lesion", rowGroups = groups(1:3); end
    rowLabel = 'Tract-based';
    if row == 1, rowLabel = 'Classical'; end
    text(ax,20,cy,{rowLabel,char(tissue)},'Rotation',90, ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    for mi = 1:numel(metrics)
        cx = centersX(mi);
        lim = radialLimits.(rowKeys(row)).(metrics(mi));
        span = diff(lim);
        for fraction = .25:.25:.75
            plot(ax,cx+radius*fraction*cos(circle),cy+radius*fraction*sin(circle), ...
                ':','Color',[.65 .65 .65]);
        end
        plot(ax,cx+radius*cos(circle),cy+radius*sin(circle),'k-');
        for ci = 1:4
            plot(ax,cx+[0 radius*cos(theta(ci))],cy+[0 radius*sin(theta(ci))], ...
                ':','Color',[.65 .65 .65]);
        end
        for gi = 1:numel(rowGroups)
            means = nan(1,4); sems = nan(1,4);
            for ci = 1:4
                ix = D.Framework == framework & D.Tissue == tissue & D.Metric == metrics(mi) & ...
                    D.Region == regions{row}(ci) & D.Subtype == rowGroups(gi);
                assert(nnz(ix) == 1,'Missing or duplicate aggregate cell.');
                means(ci) = D.Mean(ix); sems(ci) = D.SEM(ix);
            end
            center = radius*([means means(1)]-lim(1))/span;
            spread = radius*[sems sems(1)]/span;
            upper = center+spread; lower = center-spread;
            patch(ax,cx+[upper.*cos(theta) fliplr(lower.*cos(theta))], ...
                cy+[upper.*sin(theta) fliplr(lower.*sin(theta))],colors(gi,:), ...
                'FaceAlpha',.2,'EdgeColor','none');
            plot(ax,cx+center.*cos(theta),cy+center.*sin(theta),'Color',colors(gi,:));
        end
        for ci = 1:4
            labelY = 39*sign(sin(theta(ci)));
            if row == 1 && ci <= 2, labelY = 43; end
            text(ax,cx+39*sign(cos(theta(ci))),cy+labelY,labels{row}{ci}, ...
                'HorizontalAlignment','center','BackgroundColor','w','Margin',.1);
        end
        for fraction = .25:.25:.75
            text(ax,cx,cy+radius*fraction,sprintf('%.2f',lim(1)+fraction*span), ...
                'HorizontalAlignment','center','BackgroundColor','w','Margin',.1);
        end
        if row == numel(tissues)
            text(ax,cx,cy-57,metrics(mi)+" (unitless)",'HorizontalAlignment','center');
        end
    end
end

rectangle(ax,'Position',[20 12 218 17],'EdgeColor','k');
keyOrder = [1 3 2 4];
for entry = 1:numel(keyOrder)
    gi = keyOrder(entry);
    x = 24+(entry-1)*53;
    plot(ax,[x x+12],[20.5 20.5],'Color',colors(gi,:));
    text(ax,x+15,20.5,groups(gi),'Interpreter','none');
end
print(fig,outputPdf,'-dpdf','-painters');
if ~isempty(limitsJson)
    fid = fopen(limitsJson,'w');
    assert(fid >= 0,'Cannot write radial-limits JSON.');
    closeFile = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(radialLimits,'PrettyPrint',true));
end
end
