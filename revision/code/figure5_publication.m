function figure5_publication(aggregateFile,outputDir)
% Render the frozen original Figure 5 means +/- SEM, without fitting models.
% aggregateFile is the externally audited figure5_plotted_values_mounted.csv.
% It contains aggregate-only source summaries; input hashes and subjects were
% verified against the original GroupTract and icometrix workbooks before use.
% Tract workbook bytes differ from the mounted audit, but all 133 subject IDs,
% validity masks and 3,432 plotted subject-feature values agree exactly.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'helpers'));
if nargin<2,outputDir=fullfile(root,'revision','results','figures');end
D=readtable(aggregateFile,'TextType','string');
assert(height(D)==160,'Expected 160 original Figure 5 aggregate cells.');
colors=fig5_palette(); types={'RRMS','SPMS','PPMS','control'};
metrics={'T1','MTR','FA','MD'};
regions={{'periventricular','juxtacortical','infratentorial','deepwhitematter'}, ...
    {'Association','Cerebellar','Occipitoparietal','ProjectionBrainstem'}};
labels={{'PV','JC','IT','DW'},{'A','C','OP','BP'}};
% Original displayed radial limits; MD is shown in 10^-3 mm^2 s^-1.
limits={[.8 1.8],[20 44],[.05 .4],[.7 1.75]; ...
    [.8 1.45],[30 46],[.15 .48],[.7 1.4]};
units={'T1 (seconds)','MTR (%)','FA (unitless)','MD (10^{-3} mm^{2} s^{-1})'};
rowTitles={'Classical lesion','Tract-based lesion','Tract-based NAWM'};
sizeIn=[169.8 209.8]/25.4;
fig=figure('Visible','off','Color','w','Units','inches','Position',[.5 .5 sizeIn], ...
    'PaperUnits','inches','PaperSize',sizeIn,'PaperPosition',[0 0 sizeIn]);
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
annotation(fig,'textbox',[.04 .95 .94 .04],'String', ...
    'Classical versus tract-based composition across MS subtypes', ...
    'EdgeColor','none','FontName','Arial','FontSize',10,'FontWeight','bold');
for row=1:3
    k=1+(row>1);framework='Tract';if row==1,framework='Classical';end
    tissue='Lesion';if row==3,tissue='NAWM';end
    annotation(fig,'textbox',[.035 .915-(row-1)*.258 .93 .025], ...
        'String',rowTitles{row},'EdgeColor','none','FontName','Arial', ...
        'FontSize',8.5,'FontWeight','bold');
    for mi=1:4
        ax=axes(fig,'Position',[.045+(mi-1)*.241 .72-(row-1)*.258 .22 .18]);
        hold(ax,'on');axis(ax,'equal');axis(ax,'off');
        xlim(ax,[-1.28 1.28]);ylim(ax,[-1.23 1.23]);
        lim=limits{k,mi};span=diff(lim);scale=1;if mi==4,scale=1000;end
        theta=(0:4)*pi/2+pi/4;t=linspace(0,2*pi,180);
        for ring=1:5
            rr=ring/5;plot(ax,rr*cos(t),rr*sin(t),'--','Color',[.65 .65 .65],'LineWidth',.35);
        end
        plot(ax,cos(t),sin(t),'k-','LineWidth',.65);
        for ri=1:4
            plot(ax,[0 cos(theta(ri))],[0 sin(theta(ri))],'--','Color',[.65 .65 .65],'LineWidth',.35);
            text(ax,1.14*cos(theta(ri)),1.14*sin(theta(ri)),labels{k}{ri}, ...
                'FontName','Arial','FontSize',8,'FontWeight','bold', ...
                'HorizontalAlignment','center','BackgroundColor','w','Margin',.1);
        end
        nt=3+(row==3);
        for ti=1:nt
            sourceType=types{ti};if ti==4,sourceType='Control';end
            means=zeros(1,4);se=means;
            for ri=1:4
                ix=D.Framework==framework & D.Tissue==tissue & D.Metric==metrics{mi} ...
                    & D.Region==regions{k}{ri} & D.Subtype==sourceType;
                assert(sum(ix)==1,'Missing or duplicate original Figure 5 cell.');
                means(ri)=D.Mean(ix)*scale;se(ri)=D.SEM(ix)*scale;
            end
            m=([means means(1)]-lim(1))/span;
            hi=([means+se means(1)+se(1)]-lim(1))/span;
            lo=([means-se means(1)-se(1)]-lim(1))/span;
            c=colors.(types{ti})/255;
            patch(ax,[hi.*cos(theta) fliplr(lo.*cos(theta))], ...
                [hi.*sin(theta) fliplr(lo.*sin(theta))],c,'FaceAlpha',.20,'EdgeColor','none');
            plot(ax,m.*cos(theta),m.*sin(theta),'Color',c,'LineWidth',1.0);
        end
        for ring=1:4
            v=lim(1)+span*ring/5;
            if mi==2,label=sprintf('%.1f',v);elseif mi==1&&k==1,label=sprintf('%.1f',v);else,label=sprintf('%.2f',v);end
            text(ax,0,ring/5,label,'FontName','Arial','FontSize',7.2, ...
                'HorizontalAlignment','center','BackgroundColor','w','Margin',.1);
        end
        if row==3
            text(ax,0,-1.37,units{mi},'FontName','Arial','FontSize',9, ...
                'HorizontalAlignment','center','Interpreter','tex');
        end
    end
end
% Shared key is below the original 3 x 4 arrangement, at readable size.
ax=axes(fig,'Position',[.1 .115 .8 .025]);hold(ax,'on');axis(ax,'off');
h=gobjects(1,4);for ti=1:4,h(ti)=plot(ax,nan,nan,'Color',colors.(types{ti})/255,'LineWidth',1.5);end
legend(ax,h,{'RRMS','SPMS','PPMS','Control'},'Orientation','horizontal', ...
    'Location','north','Box','off','FontName','Arial','FontSize',8);
annotation(fig,'textbox',[.045 .025 .94 .075],'String', ...
    {'PV: Periventricular; JC: Juxtacortical; IT: Infratentorial; DW: Deep white matter', ...
    'A: Association; C: Cerebellar; OP: Occipitoparietal; BP: Brainstem/Projection', ...
    'Lines: mean; shaded bands: SEM. Controls are included in NAWM panels only.'}, ...
    'FontName','Arial','FontSize',7.2,'EdgeColor','none','Interpreter','none');
print(fig,fullfile(outputDir,'Figure5.pdf'),'-dpdf','-painters');
end
