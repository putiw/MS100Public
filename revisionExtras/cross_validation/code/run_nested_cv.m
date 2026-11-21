function results = run_nested_cv(varargin)
% RUN_NESTED_CV  Nested cross-validation of tract-based vs classical ridge models.
%
% Purpose
%   The manuscript reports apparent (in-sample) AUC and R2 for models with
%   16 (tract-based) and 4 (classical) predictors. This script estimates
%   out-of-sample performance for both frameworks, and for a 4-predictor
%   tract-based model, so that the between-framework difference can be
%   separated from model flexibility. It does not change any reported
%   in-sample result.
%
% Usage
%   run_nested_cv('statsDir','/path/to/derivatives/stats')
%   run_nested_cv('statsDir',d,'metrics',{'T1','MTR'},'reps',5)   % quick test
%
% Inputs read from statsDir (participant-level; never commit these)
%   clinicalScore.xlsx
%   GroupTract<metric>_All.xlsx        for T1, MTR, FA, MD
%   icometrixLesionStats.xlsx          classical lesion means for T1, MTR, FA, MD
%   GroupTract<metric>_All.csv         for NDI, ODI, FWF (FWF = ISOVF)
%   ClassicalRegionNODDI_All.csv       classical lesion means for NDI, ODI, FWF
%
% Output (aggregate only; safe to commit)
%   <outputDir>/nested_cv_results_matlab.csv
%
% Design (must stay identical to code/python_reference for the comparison)
%   Cohort      MS patients only; joint complete cases for the tract-based and
%               classical predictors of the metric (NaN and 0 = missing).
%   Outcomes    EDSS >= 4.0, MSPro >= 2 (ridge logistic);
%               log T25FW, log 9HPT-D, log 9HPT-ND, SDMT correct, MSFC-SDMT (ridge linear).
%   Models      tract16      4 clusters x L/R x (Tail, NAWM)      [manuscript model]
%               classical4   4 compartment lesion means           [manuscript model]
%               tract4_tail  4 bilateral cluster tail percentiles [dimension-matched]
%               tract8_tail  4 clusters x L/R tail percentiles
%   Outer loop  repeated K-fold (default 30 x 10), stratified for binary outcomes;
%               the same outer folds are used for every model of an outcome.
%   Inner loop  K-fold (default 5) on the training part only, selecting lambda
%               by held-out DEVIANCE (binary) or squared error (continuous).
%               Predictors are z-scored with training-fold mean and SD.
%   Measures    cv AUC  = mean over repeats of the AUC of pooled out-of-fold scores
%               cv R2   = mean over repeats of 1 - SSE_oof / SSE_(training-fold mean)
%               apparent = same model fitted to all data (lambda by 10-fold deviance)
%   Paired difference  out-of-fold scores averaged over repeats per participant,
%               then a participant bootstrap (class-stratified for binary outcomes)
%               of the difference; 2.5th-97.5th percentile interval.
%
% Requires Statistics and Machine Learning Toolbox (fitclinear, fitrlinear, cvpartition).

p = inputParser;
addParameter(p,'statsDir','',@(x) ischar(x) || isstring(x));
addParameter(p,'outputDir',fullfile(fileparts(fileparts(mfilename('fullpath'))),'results'),@(x) ischar(x) || isstring(x));
addParameter(p,'metrics',{'T1','MTR','FA','MD','NDI','ODI','FWF'},@iscell);
addParameter(p,'reps',30,@isnumeric);
addParameter(p,'outerK',10,@isnumeric);
addParameter(p,'innerK',5,@isnumeric);
addParameter(p,'lambdaGrid',logspace(-4,4,25),@isnumeric);
addParameter(p,'nBoot',2000,@isnumeric);
addParameter(p,'seed',42,@isnumeric);
parse(p,varargin{:});
o = p.Results;
statsDir = char(o.statsDir);
assert(~isempty(statsDir) && isfolder(statsDir), ...
    'Set ''statsDir'' to the folder holding clinicalScore.xlsx and the GroupTract tables.');
if ~isfolder(o.outputDir), mkdir(o.outputDir); end

clusters = {'Association','Cerebellar','Occipitoparietal','ProjectionBrainstem'};
regions  = {'periventricular','juxtacortical','infratentorial','deepwhitematter'};

% label, column, threshold, isBinary, logTransform
outcomes = { ...
    'EDSS>=4',   'EDSS',        4.0, true,  false; ...
    'MSPro>=2',  'MSPro',       2.0, true,  false; ...
    'T25FW',     'T25FW',       NaN, false, true;  ...
    '9HPT-D',    'x9HPTD',      NaN, false, true;  ...
    '9HPT-ND',   'x9HPTND',     NaN, false, true;  ...
    'SDMT',      'SDMTcorrect', NaN, false, false; ...
    'MSFC-SDMT', 'MSFC',        NaN, false, false};
modelNames = {'tract16','classical4','tract4_tail','tract8_tail'};
pairs = {'tract16','classical4'; 'tract4_tail','classical4'; 'tract8_tail','classical4'};

rows = {};
for mi = 1:numel(o.metrics)
    metric = o.metrics{mi};
    T = load_metric_table(statsDir, metric);
    T = T(~startsWith(T.SubjectID,'sub-C'), :);          % patients only

    P = struct();
    P.tract16 = {}; P.tract8_tail = {};
    for ci = 1:numel(clusters)
        for kind = {'Tail','NAWM'}
            for hemi = {'L','R'}
                P.tract16{end+1} = [clusters{ci} hemi{1} '_' kind{1}]; %#ok<AGROW>
            end
        end
        for hemi = {'L','R'}
            P.tract8_tail{end+1} = [clusters{ci} hemi{1} '_Tail']; %#ok<AGROW>
        end
    end
    P.classical4  = cellfun(@(r) [r metric], regions, 'UniformOutput', false);
    P.tract4_tail = cellfun(@(c) [c '_Tail'], clusters, 'UniformOutput', false);
    cohortCols = [P.tract16, P.classical4, P.tract4_tail];

    for oi = 1:size(outcomes,1)
        [label, col, thr, isBin, doLog] = outcomes{oi,:};
        yraw = double(T.(col));
        if doLog, yraw(yraw <= 0) = NaN; yraw = log(yraw); end
        M = T{:, cohortCols};
        keep = isfinite(yraw) & all(isfinite(M),2) & all(M ~= 0,2);
        S = T(keep,:);
        if isBin, y = double(yraw(keep) >= thr); else, y = yraw(keep); end
        n = numel(y);

        oofMean = struct();
        for ki = 1:numel(modelNames)
            name = modelNames{ki};
            X = S{:, P.(name)};
            tStart = tic;
            % Same seed before every model => identical outer folds across models.
            rng(o.seed,'twister');
            lamApp = select_lambda(X, y, isBin, 10, o.lambdaGrid);
            [b0, b] = fit_ridge(zscore(X), y, isBin, lamApp);
            lp = b0 + zscore(X) * b;
            rng(o.seed,'twister');
            [oof, base] = nested_cv(X, y, isBin, o.reps, o.outerK, o.innerK, o.lambdaGrid);
            per = zeros(o.reps,1);
            for r = 1:o.reps
                if isBin
                    per(r) = auc_rank(y, oof(r,:)');
                else
                    per(r) = 1 - sum((y - oof(r,:)').^2) / sum((y - base(r,:)').^2);
                end
            end
            if isBin
                app = auc_rank(y, lp); measure = 'AUC'; npos = sum(y);
            else
                c = corrcoef(lp, y); app = c(1,2)^2; measure = 'R2'; npos = NaN;
            end
            oofMean.(name) = mean(oof,1)';
            rows(end+1,:) = {metric, label, name, size(X,2), n, npos, measure, ...
                app, mean(per), std(per,1), app - mean(per), NaN, NaN}; %#ok<AGROW>
            fprintf('%-4s %-10s %-12s N=%3d  apparent=%.3f  cv=%.3f  (%.0fs)\n', ...
                metric, label, name, n, app, mean(per), toc(tStart));
        end

        for qi = 1:size(pairs,1)
            a = oofMean.(pairs{qi,1}); bb = oofMean.(pairs{qi,2});
            rng(o.seed + 1,'twister');
            [pe, ci] = paired_difference(y, a, bb, isBin, o.nBoot);
            if isBin, measure = 'cvAUC diff'; else, measure = 'cvR2 diff'; end
            rows(end+1,:) = {metric, label, sprintf('DELTA %s - %s', pairs{qi,1}, pairs{qi,2}), ...
                NaN, n, NaN, measure, NaN, pe, NaN, NaN, ci(1), ci(2)}; %#ok<AGROW>
        end
    end
end

results = cell2table(rows, 'VariableNames', {'metric','outcome','model','n_pred','N','pos', ...
    'measure','apparent','cv','cv_sd_across_repeats','optimism','ci_lo','ci_hi'});
outFile = fullfile(o.outputDir,'nested_cv_results_matlab.csv');
writetable(results, outFile);
fprintf('\nWrote %s (%d rows). MATLAB %s.\n', outFile, height(results), version);
end

% ------------------------------------------------------------------ data
function T = load_metric_table(statsDir, metric)
clin = readtable(fullfile(statsDir,'clinicalScore.xlsx'));
if ismember(metric, {'NDI','ODI','FWF'})
    tract = readtable(fullfile(statsDir, ['GroupTract' metric '_All.csv']));
    cls   = readtable(fullfile(statsDir, 'ClassicalRegionNODDI_All.csv'));
else
    tract = readtable(fullfile(statsDir, ['GroupTract' metric '_All.xlsx']));
    cls   = readtable(fullfile(statsDir, 'icometrixLesionStats.xlsx'));
end
tract = tract(~cellfun(@isempty, tract.SubjectID), :);
[~, iu] = unique(tract.SubjectID, 'stable'); tract = tract(iu,:);
T = outerjoin(clin, tract, 'Keys','SubjectID', 'MergeKeys',true, 'Type','left');
T = outerjoin(T,    cls,   'Keys','SubjectID', 'MergeKeys',true, 'Type','left');
% MSFC-SDMT: components z-scored against the MS cohort, then averaged.
isMS = ~startsWith(T.SubjectID,'sub-C');
arm = (1 ./ T.x9HPTD + 1 ./ T.x9HPTND) / 2;
leg = -T.T25FW;
cog = double(T.SDMTcorrect);
z = @(v) (v - mean(v(isMS),'omitnan')) ./ std(v(isMS),'omitnan');
T.MSFC = (z(arm) + z(leg) + z(cog)) / 3;
end

% ------------------------------------------------------------ modelling
function [oof, base] = nested_cv(X, y, isBin, reps, outerK, innerK, grid)
n = numel(y); oof = zeros(reps,n); base = zeros(reps,n);
for r = 1:reps
    cvo = make_partition(y, outerK, isBin);
    for k = 1:outerK
        tr = training(cvo,k); te = test(cvo,k);
        lam = select_lambda(X(tr,:), y(tr), isBin, innerK, grid);
        [Xtr, Xte] = standardize(X(tr,:), X(te,:));
        [b0, b] = fit_ridge(Xtr, y(tr), isBin, lam);
        oof(r,te)  = (b0 + Xte * b)';
        base(r,te) = mean(y(tr));
    end
end
end

function lam = select_lambda(X, y, isBin, k, grid)
% Lambda minimising held-out deviance (binary) or squared error (continuous).
cvi = make_partition(y, k, isBin);
loss = zeros(1, numel(grid));
for j = 1:k
    tr = training(cvi,j); te = test(cvi,j);
    [Xtr, Xte] = standardize(X(tr,:), X(te,:));
    [b0, B] = fit_ridge(Xtr, y(tr), isBin, grid);       % B: p x numel(grid)
    eta = Xte * B + b0;                                 % n_te x numel(grid)
    if isBin
        eta = min(max(eta,-30),30);
        loss = loss + sum(log1p(exp(eta)) - y(te) .* eta, 1);
    else
        loss = loss + sum((y(te) - eta).^2, 1);
    end
end
[~, best] = min(loss);
lam = grid(best);
end

function [b0, B] = fit_ridge(Xz, y, isBin, lambda)
% Objective: mean loss + (lambda/2)*||beta||^2, intercept unpenalised.
if isBin
    mdl = fitclinear(Xz, y, 'Learner','logistic', 'Regularization','ridge', ...
        'Lambda',lambda, 'Solver','bfgs', 'ClassNames',[0 1], ...
        'BetaTolerance',1e-8, 'GradientTolerance',1e-8);
else
    mdl = fitrlinear(Xz, y, 'Learner','leastsquares', 'Regularization','ridge', ...
        'Lambda',lambda, 'Solver','bfgs', ...
        'BetaTolerance',1e-8, 'GradientTolerance',1e-8);
end
B = mdl.Beta; b0 = mdl.Bias;                            % b0 is 1 x numel(lambda)
end

function cvp = make_partition(y, k, isBin)
if isBin
    cvp = cvpartition(y, 'KFold', k, 'Stratify', true);
else
    cvp = cvpartition(numel(y), 'KFold', k);
end
end

function [Xtr, Xte] = standardize(Xtr, Xte)
mu = mean(Xtr,1); sd = std(Xtr,0,1); sd(sd == 0) = 1;
Xtr = (Xtr - mu) ./ sd; Xte = (Xte - mu) ./ sd;
end

% -------------------------------------------------------------- metrics
function a = auc_rank(y, s)
y = logical(y(:)); r = tiedrank(s(:));
n1 = sum(y); n0 = sum(~y);
a = (sum(r(y)) - n1*(n1+1)/2) / (n1*n0);
end

function [pe, ci] = paired_difference(y, a, b, isBin, nBoot)
n = numel(y); d = zeros(nBoot,1);
if isBin
    pe = auc_rank(y,a) - auc_rank(y,b);
    pos = find(y == 1); neg = find(y == 0);
    for i = 1:nBoot
        idx = [pos(randi(numel(pos),numel(pos),1)); neg(randi(numel(neg),numel(neg),1))];
        d(i) = auc_rank(y(idx),a(idx)) - auc_rank(y(idx),b(idx));
    end
else
    sst = sum((y - mean(y)).^2);
    pe = (sum((y-b).^2) - sum((y-a).^2)) / sst;
    for i = 1:nBoot
        idx = randi(n,n,1); yy = y(idx);
        d(i) = (sum((yy-b(idx)).^2) - sum((yy-a(idx)).^2)) / sum((yy - mean(yy)).^2);
    end
end
ci = quantile(d, [0.025 0.975]);
end
