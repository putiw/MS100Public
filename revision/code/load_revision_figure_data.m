function D = load_revision_figure_data(whichFigure,varargin)
% Read frozen aggregate results only; no styling, fitting or participant data.
root=fileparts(fileparts(mfilename('fullpath')));
dataDir=fullfile(root,'results','model_results');
switch upper(char(whichFigure))
    case 'FIGURE5'
        path=fullfile(dataDir,'figure5_plotted_values.csv');
        if ~isempty(varargin) && ~isempty(varargin{1}),path=varargin{1};end
        D=readtable(path,'TextType','string');
        assert(height(D)==160,'Expected 160 frozen Figure 5 cells.');
    case 'FIGURE6'
        fullFile=fullfile(dataDir,'figure6_full_models.csv');
        matchedFile=fullfile(dataDir,'figure6_matched_models.csv');
        if ~isempty(varargin) && ~isempty(varargin{1}),fullFile=varargin{1};end
        if numel(varargin)>1 && ~isempty(varargin{2}),matchedFile=varargin{2};end
        full=readtable(fullFile,'TextType','string');
        matched=readtable(matchedFile,'TextType','string');
        assert(height(full)==70 && height(matched)==70,'Expected 70 rows per model table.');
        full=parse_table(full);matched=parse_table(matched);
        D.metrics={'T1','MTR','FA','MD','NDI','ODI','ISOVF'};
        D.labels={'T25FW','9HPT-D','9HPT-ND','SDMT','MSFC-SDMT'};
        D.keys={'T25FW','x9HPTD','x9HPTND','SDMTcorrect','MSFC_SDMT'};
        for ci=1:5
            D.panels(ci).unadjusted=collect_panel_values(full,matched,D.keys{ci},D.metrics,'R2');
            D.panels(ci).adjusted=collect_panel_values(full,matched,D.keys{ci},D.metrics,'R2_Demo');
        end
    case 'S1'
        path=fullfile(dataDir,'supplementary_figure_s1_correlations.json');
        if ~isempty(varargin) && ~isempty(varargin{1}),path=varargin{1};end
        payload=jsondecode(fileread(path));
        assert(numel(payload.cells)==182,'Expected 182 correlation cells.');
        assert(payload.multiple_testing.global_family.multiplier==182,'Expected 182-test family.');
        D=build_panels(payload.cells);
    otherwise
        error('Unknown aggregate figure: %s',whichFigure);
end
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
