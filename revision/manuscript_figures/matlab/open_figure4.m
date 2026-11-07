function figures = open_figure4()
% Open the quantitative parts of Figure 4 from aggregate tables only.
repoDir=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
inputDir=fullfile(repoDir,'revisionExtras','followup_20261013','results');
barData=readtable(fullfile(inputDir,'figure4a_bar_by_subtype_caption_matched_89.csv'),'TextType','string');
radarData=readtable(fullfile(inputDir,'lesion_load_by_region_subtype_caption_matched_89.csv'),'TextType','string');
groups=["RRMS","PPMS","SPMS"];
counts=[49 17 23];
frameworks=["Classical","Tract"];
barFrameworks=["WholeBrain","Tract"];
metrics=["LN","LV","Lnorm"];
metricNames=["Lesion number (count)","Lesion volume (mL)","Normalized lesion volume (%)"];
regions={["periventricular","juxtacortical","infratentorial","deepwhitematter"], ...
         ["Association","Cerebellar","Occipitoparietal","ProjectionBrainstem"]};
regionLabels={["Periventricular","Juxtacortical","Infratentorial","Deep WM"], ...
              ["Association","Cerebellar","Occipitoparietal","Projection/Brainstem"]};

assert(height(barData)==6 && height(radarData)==72,'Unexpected Figure 4 table size.');
means=zeros(3,2);
sds=zeros(3,2);
for g=1:3
    for f=1:2
        row=barData(barData.Subtype==groups(g) & barData.Framework==barFrameworks(f),:);
        assert(height(row)==1 && row.N==counts(g),'Unexpected Figure 4A aggregate.');
        means(g,f)=row.Mean;
        sds(g,f)=row.SD;
    end
end
figures(1)=figure('Name','Figure 4A');
ax=axes(figures(1));hold(ax,'on');
bars=bar(ax,means);
drawnow;
for f=1:2
    errorbar(ax,bars(f).XEndPoints,means(:,f),sds(:,f),'k','LineStyle','none','HandleVisibility','off');
end
xticks(ax,1:3);
xticklabels(ax,compose('%s (n=%d)',groups(:),counts(:)));
ylabel(ax,'Lesion volume (mL)');
ylim(ax,[0 50]);
legend(ax,bars,{'Whole brain','Sum of four tract regions'},'Location','northwest');
title(ax,'Mean ± SD');

figures(2)=figure('Name','Figure 4B');
layout=tiledlayout(figures(2),2,3);
theta=(0:4)*pi/2+pi/4;
for f=1:2
    for m=1:3
        ax=polaraxes(layout);
        ax.Layout.Tile=(f-1)*3+m;
        hold(ax,'on');
        mu=zeros(4,3);se=mu;
        for g=1:3
            for k=1:4
                row=radarData(radarData.Framework==frameworks(f) & radarData.Metric==metrics(m) ...
                    & radarData.Region==regions{f}(k) & radarData.Subtype==groups(g),:);
                assert(height(row)==1 && row.N==counts(g),'Unexpected Figure 4B aggregate.');
                mu(k,g)=row.Mean;
                se(k,g)=row.SEM;
            end
        end
        if metrics(m)=="LV"
            mu=mu/1000;
            se=se/1000;
        end
        colors=ax.ColorOrder;
        lines=gobjects(1,3);
        for g=1:3
            center=mu([1:4 1],g);
            spread=se([1:4 1],g);
            lines(g)=polarplot(ax,theta,center,'Color',colors(g,:));
            polarplot(ax,theta,center+spread,':','Color',colors(g,:),'HandleVisibility','off');
            polarplot(ax,theta,max(0,center-spread),':','Color',colors(g,:),'HandleVisibility','off');
        end
        thetaticks(ax,45:90:315);thetaticklabels(ax,regionLabels{f});
        rlim(ax,[0 1.1*max(mu+se,[],'all')]);
        title(ax,frameworks(f)+" "+metricNames(m));
        if f==1 && m==1
            legend(ax,lines,groups,'Location','southoutside','Orientation','horizontal');
        end
    end
end
figures(1).UserData=table2struct(barData);
figures(2).UserData=table2struct(radarData);
end
