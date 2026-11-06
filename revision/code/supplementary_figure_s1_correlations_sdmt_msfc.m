function [report,fig] = supplementary_figure_s1_correlations_sdmt_msfc(varargin)
% Plain MATLAB S1 from frozen aggregate correlations, without model fitting.
p=inputParser;root=fileparts(fileparts(mfilename('fullpath')));
addParameter(p,'resultsFile',fullfile(root,'results','model_results','supplementary_figure_s1_correlations.json'));
addParameter(p,'outputDir',fullfile(root,'results','figures'));
addParameter(p,'previewDir',fullfile(root,'work','figure_previews'));
addParameter(p,'doExport',true);
addParameter(p,'keepOpen',false);parse(p,varargin{:});
D=load_revision_figure_data('S1',p.Results.resultsFile);
pdfPath=fullfile(p.Results.outputDir,'Supplementary_FigureS1.pdf');
if p.Results.doExport
    if ~isfolder(p.Results.outputDir),mkdir(p.Results.outputDir);end
    if ~isfolder(p.Results.previewDir),mkdir(p.Results.previewDir);end
end
fig=gobjects(0);
if p.Results.doExport || p.Results.keepOpen
    sz=[169.8 209.8]/25.4;fig=figure('Visible','off','Units','inches', ...
        'Position',[.5 .5 sz],'PaperUnits','inches','PaperSize',sz,'PaperPosition',[0 0 sz]);
    tl=tiledlayout(fig,3,2);palette=parula(256);
    for pi=1:6
        ax=nexttile(tl);hold(ax,'on');A=D(pi).corr_matrix;
        for i=1:4
            for j=1:7
                v=A(i,j);c=palette(min(256,max(1,floor(v*256)+1)),:);
                patch(ax,j+[-.5 .5 .5 -.5],i+[-.5 -.5 .5 .5],c,'EdgeColor','none');
                label=sprintf('%.2f',v);if D(pi).significant(i,j),label=[label '*'];end
                text(ax,j,i,label,'HorizontalAlignment','center','FontSize',8);
            end
        end
        set(ax,'XLim',[.5 7.5],'YLim',[.5 4.5],'YDir','reverse','XTick',1:7, ...
            'XTickLabel',D(pi).clinical_labels,'YTick',1:4,'YTickLabel',region_labels(pi),'FontSize',8);
        xtickangle(ax,55);clim(ax,[0 1]);colormap(ax,palette);colorbar(ax);
        title(ax,D(pi).title,'FontSize',8,'Interpreter','none');
    end
    if p.Results.doExport
        print(fig,pdfPath,'-dpdf','-painters');
        print(fig,fullfile(p.Results.previewDir,'Supplementary_FigureS1.png'),'-dpng','-r300');
    end
    if p.Results.keepOpen
        fig.Name='Supplementary Figure S1';fig.UserData=D;fig.Visible='on';
    else
        close(fig);fig=gobjects(0);
    end
end
report=struct('figure',pdfPath,'panels',6,'correlation_cells',182, ...
    'bonferroni_multiplier',182,'exported',p.Results.doExport);
report.source_script='code/supplementary_figure_s1_correlations_sdmt_msfc.m';
report.results_file=char(p.Results.resultsFile);
end

function labels=region_labels(panel)
if mod(panel,2)==1,labels={'A','C','OP','BP'};else,labels={'PV','JC','IT','DW'};end
end
