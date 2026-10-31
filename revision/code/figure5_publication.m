function figure5_publication(aggregateFile,outputDir)
% Plain MATLAB Figure 5 from tracked aggregate inputs; no manuscript styling.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,aggregateFile=[];end
if nargin<2 || isempty(outputDir),outputDir=fullfile(root,'results','figures');end
if ~isfolder(outputDir),mkdir(outputDir);end
D=load_revision_figure_data('Figure5',aggregateFile);
metrics={'T1','MTR','FA','MD'};types={'RRMS','SPMS','PPMS','Control'};
regions={{'periventricular','juxtacortical','infratentorial','deepwhitematter'}, ...
    {'Association','Cerebellar','Occipitoparietal','ProjectionBrainstem'}};
labels={{'PV','JC','IT','DW'},{'A','C','OP','BP'}};
limits={[.8 1.8],[20 44],[.05 .4],[.7 1.75]; ...
    [.8 1.45],[30 46],[.15 .48],[.7 1.4]};
units={'T1 (seconds)','MTR (%)','FA (unitless)','MD (10^{-3} mm^{2} s^{-1})'};
sizeIn=[169.8 209.8]/25.4;
fig=figure('Visible','off','Units','inches','Position',[.5 .5 sizeIn], ...
    'PaperUnits','inches','PaperSize',sizeIn,'PaperPosition',[0 0 sizeIn]);
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
for row=1:3
    k=1+(row>1);fw='Tract';if row==1,fw='Classical';end
    tissue='Lesion';if row==3,tissue='NAWM';end
    for mi=1:4
        ax=subplot(3,4,(row-1)*4+mi,'Parent',fig);hold(ax,'on');axis(ax,'equal');
        theta=(0:4)*pi/2+pi/4;t=linspace(0,2*pi,120);lim=limits{k,mi};span=diff(lim);
        for ring=1:5,plot(ax,ring/5*cos(t),ring/5*sin(t),':','Color',[.7 .7 .7],'HandleVisibility','off');end
        h=gobjects(1,3+(row==3));cs=ax.ColorOrder;
        for ti=1:numel(h)
            means=zeros(1,4);se=means;
            for ri=1:4
                ix=D.Framework==fw & D.Tissue==tissue & D.Metric==metrics{mi} ...
                    & D.Region==regions{k}{ri} & D.Subtype==types{ti};
                assert(sum(ix)==1,'Missing/duplicate aggregate cell.');
                scale=1;if mi==4,scale=1000;end
                means(ri)=D.Mean(ix)*scale;se(ri)=D.SEM(ix)*scale;
            end
            m=([means means(1)]-lim(1))/span;
            hi=([means+se means(1)+se(1)]-lim(1))/span;
            lo=([means-se means(1)-se(1)]-lim(1))/span;
            patch(ax,[hi.*cos(theta) fliplr(lo.*cos(theta))], ...
                [hi.*sin(theta) fliplr(lo.*sin(theta))],cs(ti,:), ...
                'FaceAlpha',.2,'EdgeColor','none','HandleVisibility','off');
            h(ti)=plot(ax,m.*cos(theta),m.*sin(theta),'Color',cs(ti,:));
        end
        for ri=1:4,text(ax,1.14*cos(theta(ri)),1.14*sin(theta(ri)),labels{k}{ri}, ...
                'HorizontalAlignment','center','FontSize',8);end
        for ring=1:4,text(ax,0,ring/5,sprintf('%.2f',lim(1)+span*ring/5), ...
                'HorizontalAlignment','center','FontSize',8);end
        axis(ax,'off');xlim(ax,[-1.35 1.35]);ylim(ax,[-1.35 1.35]);
        title(ax,{[fw ' ' tissue],metrics{mi}},'FontSize',8);
        if row==3,text(ax,0,-1.6,units{mi},'HorizontalAlignment','center','FontSize',9);end
        if row==3 && mi==1,legend(ax,h,types,'Location','southoutside','FontSize',8);end
    end
end
print(fig,fullfile(outputDir,'Figure5.pdf'),'-dpdf','-painters');
end
