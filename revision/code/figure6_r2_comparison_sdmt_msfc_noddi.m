function [report,figures] = figure6_r2_comparison_sdmt_msfc_noddi(varargin)
% Plain MATLAB Figure 6/S2, using frozen model results; no model fitting.
p=inputParser;root=fileparts(fileparts(mfilename('fullpath')));
addParameter(p,'fullFile',fullfile(root,'results','model_results','figure6_full_models.csv'));
addParameter(p,'sameNFile',fullfile(root,'results','model_results','figure6_matched_models.csv'));
addParameter(p,'outputDir',fullfile(root,'results','figures'));
addParameter(p,'previewDir',fullfile(root,'work','figure_previews'));
addParameter(p,'doExport',true);
addParameter(p,'keepOpen',false);
addParameter(p,'variants',1:2,@(x) isnumeric(x) && all(ismember(x,1:2)));
parse(p,varargin{:});
D=load_revision_figure_data('Figure6',p.Results.fullFile,p.Results.sameNFile);
names={'Figure6','Supplementary_FigureS2'};fields={'unadjusted','adjusted'};
titles={'Unadjusted R²','R² adjusted for age and recorded gender'};
if p.Results.doExport
    if ~isfolder(p.Results.outputDir),mkdir(p.Results.outputDir);end
    if ~isfolder(p.Results.previewDir),mkdir(p.Results.previewDir);end
end
figures=gobjects(0,1);
for fi=p.Results.variants(:).'
    sz=[169.8 209.8]/25.4;
    fig=figure('Visible','off','Units','inches','Position',[.5 .5 sz], ...
        'PaperUnits','inches','PaperSize',sz,'PaperPosition',[0 0 sz]);
    tl=tiledlayout(fig,5,1);title(tl,titles{fi});first=[];
    for ci=1:5
        v=D.panels(ci).(fields{fi});ax=nexttile(tl);hold(ax,'on');
        b=bar(ax,[v.r2Tract(:),v.r2Classic(:)]);drawnow;
        for mi=1:7
            vals=[v.r2Tract(mi),v.r2Classic(mi)];matched=[v.r2TractSameN(mi),v.r2ClassicSameN(mi)];ns=[v.nTract(mi),v.nClassic(mi)];
            for model=1:2
                x=b(model).XEndPoints(mi);
                plot(ax,[x-.15 x+.15],[matched(model) matched(model)],'k-','HandleVisibility','off');
                text(ax,x,max(vals(model),matched(model))+.018,sprintf('%d',ns(model)), ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
            end
        end
        xticks(ax,1:7);xticklabels(ax,D.metrics);ylim(ax,[0 .56]);xlim(ax,[.45 7.55]);
        ylabel(ax,'R²');title(ax,D.labels{ci});ax.FontSize=8;
        if ci==1
            first=ax;h=gobjects(3,1);h(1:2)=b(:);h(3)=plot(ax,nan,nan,'k-');
        end
    end
    legend(first,h,{'Tract-based','Classical','Matched MS cohort'},'Location','northoutside','Orientation','horizontal');
    xlabel(tl,'MRI metric');
    if p.Results.doExport
        print(fig,fullfile(p.Results.outputDir,[names{fi} '.pdf']),'-dpdf','-painters');
        print(fig,fullfile(p.Results.previewDir,[names{fi} '.png']),'-dpng','-r300');
    end
    if p.Results.keepOpen
        fig.Name=names{fi};fig.UserData=D;fig.Visible='on';
        figures(end+1,1)=fig; %#ok<AGROW>
    else
        close(fig);
    end
end
report=struct('figure6',fullfile(p.Results.outputDir,'Figure6.pdf'), ...
    'supplementary_figure_s2',fullfile(p.Results.outputDir,'Supplementary_FigureS2.pdf'), ...
    'layout','five rows by one column','exported',p.Results.doExport);
report.source_script='code/figure6_r2_comparison_sdmt_msfc_noddi.m';
report.metrics=D.metrics;report.outcomes=D.labels;
report.noddi_framework='tract-based and classical-region';
end
