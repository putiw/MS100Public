function fig = open_figure1(includeAddedScores)
% Open the clinical histograms and radar from frozen aggregates.
% open_figure1(true) adds SDMT and MSFC-SDMT in a separate version.
arguments
    includeAddedScores (1,1) logical = false
end
here=fileparts(mfilename('fullpath'));
D=jsondecode(fileread(fullfile(here,'figure1_aggregates.json')));
maxima=D.radar.maxima(:).';
minima=zeros(size(maxima));
if includeAddedScores
    added=jsondecode(fileread(fullfile(here,'figure1_sdmt_msfc_aggregates.json')));
    D.histograms=[D.histograms(:);added.histograms(:)];
    D.radar.labels=[D.radar.labels(:);{'SDMT';'MSFC-SDMT'}];
    D.radar.means=[D.radar.means,added.radar_means.'];
    D.radar.sems=[D.radar.sems,added.radar_sems.'];
    D.radar.counts=[D.radar.counts,added.sample_sizes(:,2:4).'];
    maxima=[maxima 80 1];
    minima=[minima 30 -1.5];
    fig=figure('Name','Figure 1 alternative: SDMT and MSFC-SDMT');
    layout=tiledlayout(fig,3,3);
else
    fig=figure('Name','Figure 1');
    layout=tiledlayout(fig,2,3);
end
for k=1:numel(D.histograms)
    h=D.histograms(k);
    ax=nexttile(layout);
    bars=bar(ax,h.bins,h.counts,'stacked');
    ylabel(ax,'Participant count');
    switch h.label
        case {'T25FW','9HPT-D','9HPT-ND'}
            xlabel(ax,'Time (seconds)');
        case 'SDMT'
            xlabel(ax,'Correct responses');
        case 'MSFC-SDMT'
            xlabel(ax,'Composite z-score');
        otherwise
            xlabel(ax,'Score');
    end
    title(ax,h.label);
    ylim(ax,[0 h.yMax]);
    if k==1
        legend(ax,bars,D.groups,'Location','best');
    end
end
ax=polaraxes(layout);
ax.Layout.Tile=numel(D.histograms)+1;
if includeAddedScores
    ax.Layout.TileSpan=[1 2];
end
hold(ax,'on');
n=numel(D.radar.labels);
theta=(0:n)*2*pi/n+pi/2;
series=gobjects(1,3);
for g=1:3
    mu=(D.radar.means(g,:)-minima)./(maxima-minima);
    if includeAddedScores
        mu(end-1:end)=1-mu(end-1:end);
    end
    se=D.radar.sems(g,:)./(maxima-minima);
    color=ax.ColorOrder(g,:);
    order=[1:n 1];
    series(g)=polarplot(ax,theta,mu(order),'Color',color);
    polarplot(ax,theta,mu(order)+se(order),':','Color',color,'HandleVisibility','off');
    polarplot(ax,theta,mu(order)-se(order),':','Color',color,'HandleVisibility','off');
end
[ticks,labelOrder]=sort(mod(rad2deg(theta(1:n)),360));
thetaticks(ax,ticks);
thetaticklabels(ax,D.radar.labels(labelOrder));
title(ax,'Mean ± SEM (normalized)');
legend(ax,series,D.groups(2:4),'Location','southoutside');
fig.UserData=D;
end
