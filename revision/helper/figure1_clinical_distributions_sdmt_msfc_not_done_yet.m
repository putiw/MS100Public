function report = figure1_clinical_distributions_sdmt_msfc(varargin)
%FIGURE1_CLINICAL_DISTRIBUTIONS_SDMT_MSFC Reproduce revised Figure 1 parts.
%
% Copied from MS100Public/fig1_clinical_distributions.m. The changes are:
%   1. add the SDMT-based MSFC score;
%   2. show seven histograms in a 4-over-3 layout; and
%   3. show seven measures in the radar plot.
%
% The function leaves two MATLAB figures open for assembly in Illustrator
% and also exports each figure as a separate vector PDF. clinicalScore.xlsx
% is read-only; MSFC-SDMT is calculated in memory using the same definition
% as generate_manuscript_tables.m and run_all_qmri_manuscript_models.m.
%
% Example:
%   figure1_clinical_distributions_sdmt_msfc( ...
%       'clinicalFile', '/path/to/clinicalScore.xlsx');


p = inputParser;
script_dir = fileparts(mfilename('fullpath'));
package_root = fileparts(script_dir);
public_root = fileparts(package_root);

addParameter(p, 'clinicalFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'outputDir', fullfile(package_root, 'results', 'figures'), ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'doExport', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'visible', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});

clinical_file = resolve_clinical_file(char(p.Results.clinicalFile), public_root);
output_dir = char(p.Results.outputDir);
if p.Results.visible, visible = 'on'; else, visible = 'off'; end

clinical = readtable(clinical_file);
clinical = add_msfc_sdmt(clinical);
group = numeric_vector(clinical.Group);

group_codes = [0 1 2 3];
group_names = {'Control','RRMS','PPMS','SPMS'};
expected_group_n = [43 49 17 23];
actual_group_n = arrayfun(@(code) sum(group == code), group_codes);
assert(isequal(actual_group_n, expected_group_n), ...
    'Expected group counts [%s], found [%s].', ...
    num2str(expected_group_n), num2str(actual_group_n));

config = jsondecode(fileread(fullfile(public_root, 'configs', 'config.json')));
colour_fields = {'control','RRMS','PPMS','SPMS'};
group_colours = zeros(4,3);
for i = 1:4
    group_colours(i,:) = config.msTypeColors.(colour_fields{i}) / 255;
end

fields = {'EDSS','T25FW','x9HPTD','x9HPTND', ...
    'SDMTcorrect','MSFC_SDMT','MSPro'};
labels = {'EDSS','T25FW','9HPT-D','9HPT-ND', ...
    'SDMT','MSFC-SDMT','MSPro'};
xlabels = {'Score','Time (seconds)','Time (seconds)','Time (seconds)', ...
    'Correct responses','Composite z-score','Score'};

histogram_size_mm = [109 72.5];
histogram_figure = draw_histograms(clinical, group, group_codes, ...
    group_names, group_colours, fields, labels, xlabels, ...
    histogram_size_mm, visible);

radar_fields = {'MSFC_SDMT','SDMTcorrect','x9HPTND','T25FW', ...
    'EDSS','x9HPTD','MSPro'};
radar_labels = {'MSFC-SDMT','SDMT','9HPT-ND','T25FW', ...
    'EDSS','9HPT-D','MSPro'};
radar_limits = [-1.5 1.0; 30 80; 0 40; 0 20; 0 7; 0 40; 0 3];
higher_is_better = [true true false false false false false];
patient_codes = [1 2 3];
patient_colours = group_colours(2:4,:);
[radar_means, radar_sem, radar_n] = group_summary( ...
    clinical, group, radar_fields, patient_codes);

radar_size_mm = [60 72.5];
radar_figure = draw_radar_figure(radar_means, radar_sem, radar_labels, ...
    radar_limits, higher_is_better, patient_colours, radar_size_mm, visible);

histogram_pdf = fullfile(output_dir, 'Figure1_histograms_SDMT_MSFC.pdf');
radar_pdf = fullfile(output_dir, 'Figure1_radar_SDMT_MSFC.pdf');
if p.Results.doExport
    if ~isfolder(output_dir), mkdir(output_dir); end
    export_pdf(histogram_figure, histogram_pdf, histogram_size_mm);
    export_pdf(radar_figure, radar_pdf, radar_size_mm);
end

report = struct();
report.source_script = 'code/figure1_clinical_distributions_sdmt_msfc.m';
report.clinical_file = clinical_file;
report.histogram_pdf = histogram_pdf;
report.radar_pdf = radar_pdf;
report.group_names = group_names;
report.group_n = actual_group_n;
report.radar_fields = radar_fields;
report.radar_group_names = {'RRMS','PPMS','SPMS'};
report.radar_means = radar_means;
report.radar_sem = radar_sem;
report.radar_n = radar_n;
report.msfc_complete_n = sum(isfinite(clinical.MSFC_SDMT));
report.exported = p.Results.doExport;

if ~p.Results.visible
    close(histogram_figure);
    close(radar_figure);
end
end


function fig = draw_histograms(T, group, group_codes, group_names, ...
        group_colours, fields, labels, xlabels, figure_size_mm, visible)
figure_size = figure_size_mm/25.4;
fig = figure('Visible',visible, 'Color','white', 'Units','inches', ...
    'Position',[0.5 0.5 figure_size], 'DefaultAxesFontName','Arial', ...
    'DefaultTextFontName','Arial');

left = 0.065;
plot_width = 0.895;
column_gap = 0.040;
column_width = (plot_width-3*column_gap)/4;
row_y = [0.555 0.105];
row_height = 0.305;
first_bars = gobjects(4,1);

for i = 1:numel(fields)
    if i <= 4
        row = 1;
        column = i;
        x = left+(column-1)*(column_width+column_gap);
    else
        row = 2;
        column = i-4;
        x = left+0.5*(column_width+column_gap) + ...
            (column-1)*(column_width+column_gap);
    end

    ax = axes('Parent',fig, ...
        'Position',[x row_y(row) column_width row_height]);
    values = numeric_vector(T.(fields{i}));
    finite_values = values(isfinite(values));
    assert(~isempty(finite_values), 'No finite values for %s.', fields{i});

    [~,edges] = histcounts(finite_values);
    stacked = zeros(numel(edges)-1, numel(group_codes));
    for g = 1:numel(group_codes)
        stacked(:,g) = histcounts(values(group == group_codes(g)), edges).';
    end

    % Match the original Figure 1: bars are anchored at the lower bin
    % edges, occupy 80% of the bin spacing, and retain a visible gap.
    bars = bar(ax, edges(1:end-1), stacked, 0.8, 'stacked', ...
        'EdgeColor','none');
    for g = 1:numel(group_codes)
        bars(g).FaceColor = group_colours(g,:);
        bars(g).FaceAlpha = 0.75;
    end
    if i == 1, first_bars = bars; end

    title(ax, labels{i}, 'FontSize',8.5, 'FontWeight','bold', ...
        'FontAngle','italic', 'Interpreter','none');
    xlabel(ax, xlabels{i}, 'FontSize',7.0, 'FontWeight','bold', ...
        'FontAngle','italic', 'Interpreter','none');
    ylabel(ax, 'Count', 'FontSize',7.0, 'FontWeight','bold', ...
        'FontAngle','italic');

    max_count = max(sum(stacked,2));
    upper = max(10,ceil(max_count/10)*10);
    ylim(ax,[0 upper]);
    yticks(ax,0:10:upper);
    set(ax, 'FontSize',6.4, 'FontAngle','italic', 'LineWidth',1.0, ...
        'TickDir','out', 'Box','off', 'Layer','top');
    apply_histogram_ticks(ax,fields{i});
end

lg = legend(first_bars, group_names, 'Orientation','horizontal', ...
    'Position',[0.30 0.915 0.40 0.055], 'Box','off', ...
    'FontName','Arial', 'FontSize',6.8, 'FontAngle','italic');
lg.AutoUpdate = 'off';
end


function apply_histogram_ticks(ax,field)
% Preserve the horizontal tick pattern visible in the Illustrator figure.
switch field
    case 'EDSS'
        ticks = 0:7;
        tick_font_size = 6.0;
    case 'T25FW'
        ticks = [3 9 15 21 27 45];
        tick_font_size = 5.5;
    case 'x9HPTD'
        ticks = 10:5:55;
        tick_font_size = 4.8;
    case 'x9HPTND'
        ticks = 20:20:80;
        tick_font_size = 6.0;
    case 'SDMTcorrect'
        ticks = 20:10:80;
        tick_font_size = 5.5;
    case 'MSFC_SDMT'
        ticks = -4:2:2;
        tick_font_size = 6.0;
    case 'MSPro'
        ticks = 1:3;
        tick_font_size = 6.0;
    otherwise
        return;
end
xticks(ax,ticks);
xtickangle(ax,0);
ax.FontSize = tick_font_size;
end


function fig = draw_radar_figure(means, sem, labels, limits, ...
        higher_is_better, colours, figure_size_mm, visible)
figure_size = figure_size_mm/25.4;
fig = figure('Visible',visible, 'Color','white', 'Units','inches', ...
    'Position',[0.5 0.5 figure_size], 'DefaultAxesFontName','Arial', ...
    'DefaultTextFontName','Arial');
ax_position = [0.035 0.035 0.930 0.860];
ax = axes('Parent',fig, 'Position',ax_position);
hold(ax,'on');
axis(ax,'equal');
axis(ax,'off');

n_measures = numel(labels);
theta = pi/2-(0:n_measures-1)*2*pi/n_measures;
closed_theta = [theta theta(1)];

for radius = 0.2:0.2:1
    plot(ax, radius*cos(closed_theta), radius*sin(closed_theta), ...
        'Color',[0.88 0.88 0.88], 'LineWidth',0.45, ...
        'HandleVisibility','off');
end
for i = 1:n_measures
    plot(ax, [0 cos(theta(i))], [0 sin(theta(i))], ...
        'Color',[0.76 0.76 0.76], 'LineWidth',0.55, ...
        'HandleVisibility','off');
end

normalised = nan(size(means));
for g = 1:size(means,2)
    radius = normalise_radar(means(:,g), limits, higher_is_better);
    radius_a = normalise_radar(means(:,g)+sem(:,g), limits, higher_is_better);
    radius_b = normalise_radar(means(:,g)-sem(:,g), limits, higher_is_better);
    outer = max(radius_a,radius_b);
    inner = min(radius_a,radius_b);
    patch(ax, ...
        [outer(:).'.*cos(theta),fliplr(inner(:).'.*cos(theta))], ...
        [outer(:).'.*sin(theta),fliplr(inner(:).'.*sin(theta))], ...
        colours(g,:), 'FaceAlpha',0.24, 'EdgeColor','none', ...
        'HandleVisibility','off');
    plot(ax, [radius(:).' radius(1)].*cos(closed_theta), ...
        [radius(:).' radius(1)].*sin(closed_theta), '-o', ...
        'Color',colours(g,:), 'MarkerFaceColor',colours(g,:), ...
        'MarkerEdgeColor',colours(g,:), 'MarkerSize',3.2, ...
        'LineWidth',1.3, 'HandleVisibility','off');
    normalised(:,g) = radius;
end

for i = 1:n_measures
    label_radius = 1.20;
    x = label_radius*cos(theta(i));
    y = label_radius*sin(theta(i));
    [horizontal,vertical] = text_alignment(x,y);
    text(ax,x,y,labels{i}, 'HorizontalAlignment',horizontal, ...
        'VerticalAlignment',vertical, 'FontName','Arial', 'FontSize',7.0, ...
        'FontWeight','bold', 'FontAngle','italic', 'Interpreter','none', ...
        'Clipping','off');

    radial_offsets = [0.060 0.020 0.080];
    tangent_offsets = [-0.080 -0.020 0.090];
    if abs(normalised(i,2)-normalised(i,3)) < 0.08
        radial_offsets(2:3) = [0.000 0.120];
        tangent_offsets(2:3) = [-0.100 0.130];
    end
    radial_unit = [cos(theta(i)) sin(theta(i))];
    tangent_unit = [-sin(theta(i)) cos(theta(i))];
    for g = 1:size(means,2)
        radius = min(1.07,max(0.05,normalised(i,g)+radial_offsets(g)));
        position = radius*radial_unit+tangent_offsets(g)*tangent_unit;
        if strcmp(labels{i},'MSFC-SDMT')
            value_text = sprintf('%.2f',means(i,g));
        else
            value_text = sprintf('%.1f',means(i,g));
        end
        text(ax,position(1),position(2),value_text, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'FontName','Arial', 'FontSize',4.6, 'FontWeight','bold', ...
            'FontAngle','italic', 'Color','k', 'Interpreter','none', ...
            'Clipping','off');
    end
end

legend_names = {'RRMS','PPMS','SPMS'};
legend_handles = gobjects(3,1);
for g = 1:3
    legend_handles(g) = plot(ax,nan,nan,'-o', 'Color',colours(g,:), ...
        'MarkerFaceColor',colours(g,:), 'MarkerSize',3.2, 'LineWidth',1.3);
end
lg = legend(ax,legend_handles,legend_names, 'Orientation','horizontal', ...
    'Position',[0.24 0.915 0.52 0.055], 'Box','off', ...
    'FontName','Arial', 'FontSize',6.2, 'FontAngle','italic');
lg.AutoUpdate = 'off';
ax.Position = ax_position;

xlim(ax,[-1.65 1.65]);
ylim(ax,[-1.48 1.35]);
hold(ax,'off');
end


function export_pdf(fig,path,figure_size_mm)
figure_size = figure_size_mm/25.4;
set(fig, 'PaperUnits','inches', 'PaperSize',figure_size, ...
    'PaperPosition',[0 0 figure_size], 'PaperPositionMode','manual', ...
    'InvertHardcopy','off');
print(fig,path,'-dpdf','-vector');
end


function clinical_file = resolve_clinical_file(requested,public_root)
if ~isempty(requested)
    assert(isfile(requested), 'Missing clinical workbook: %s', requested);
    clinical_file = requested;
    return;
end

config_path = fullfile(public_root,'configs','config.json');
assert(isfile(config_path), 'Missing public configuration: %s', config_path);
config = jsondecode(fileread(config_path));
clinical_file = fullfile(config.bidsDir,config.statsDir,'clinicalScore.xlsx');
assert(isfile(clinical_file), ['Missing clinical workbook: %s\n' ...
    'Mount the study volume or pass ''clinicalFile'' explicitly.'], clinical_file);
end


function T = add_msfc_sdmt(T)
required = {'SubjectID','Group','T25FW','x9HPTD','x9HPTND','SDMTcorrect'};
assert(isempty(setdiff(required,T.Properties.VariableNames)), ...
    'MSFC source columns are missing.');

subjects = string(T.SubjectID);
assert(numel(unique(subjects)) == height(T), ...
    'MSFC scoring requires one row per unique subject.');

components = [numeric_vector(T.T25FW),numeric_vector(T.x9HPTD), ...
    numeric_vector(T.x9HPTND),numeric_vector(T.SDMTcorrect)];
valid = all(isfinite(components),2) & all(components(:,1:3) > 0,2) & ...
    components(:,4) >= 0;
reference = valid & ~startsWith(subjects,'sub-C');
assert(height(T) == 132 && sum(valid) == 127 && sum(reference) == 84, ...
    'Unexpected MSFC-SDMT cohort.');

arm = nan(height(T),1);
arm(valid) = (1./components(valid,2)+1./components(valid,3))/2;
z_arm = nan(height(T),1);
z_leg = nan(height(T),1);
z_cog = nan(height(T),1);
z_arm(valid) = (arm(valid)-mean(arm(reference)))/std(arm(reference));
z_leg(valid) = -(components(valid,1)-mean(components(reference,1))) / ...
    std(components(reference,1));
z_cog(valid) = (components(valid,4)-mean(components(reference,4))) / ...
    std(components(reference,4));
score = mean([z_arm,z_leg,z_cog],2,'omitmissing');
score(~valid) = NaN;
T.MSFC_SDMT = score;
end


function [means,sem,n] = group_summary(T,group,fields,codes)
means = nan(numel(fields),numel(codes));
sem = nan(size(means));
n = zeros(size(means));
for i = 1:numel(fields)
    values = numeric_vector(T.(fields{i}));
    for g = 1:numel(codes)
        keep = group == codes(g) & isfinite(values);
        group_values = values(keep);
        n(i,g) = numel(group_values);
        assert(n(i,g) > 0, 'No data for %s, group %d.', fields{i}, codes(g));
        means(i,g) = mean(group_values);
        sem(i,g) = std(group_values)/sqrt(n(i,g));
    end
end
end


function radius = normalise_radar(values,limits,higher_is_better)
radius = (values-limits(:,1))./(limits(:,2)-limits(:,1));
radius(higher_is_better(:)) = 1-radius(higher_is_better(:));
radius = min(1,max(0,radius));
end


function [horizontal,vertical] = text_alignment(x,y)
if abs(y) < 0.35
    horizontal = 'center';
elseif x > 0.10
    horizontal = 'left';
elseif x < -0.10
    horizontal = 'right';
else
    horizontal = 'center';
end
if y > 0.10
    vertical = 'bottom';
elseif y < -0.10
    vertical = 'top';
else
    vertical = 'middle';
end
end


function values = numeric_vector(values)
if isnumeric(values) || islogical(values)
    values = double(values);
else
    values = str2double(string(values));
end
values = values(:);
end
