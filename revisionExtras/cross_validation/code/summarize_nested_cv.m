function summarize_nested_cv(resultsFile, referenceFile)
% SUMMARIZE_NESTED_CV  Headline counts and ranges from run_nested_cv output,
% and (optionally) a comparison with the Python reference results.
%
%   summarize_nested_cv()                              % MATLAB results only
%   summarize_nested_cv(resultsFile, referenceFile)    % also compare
here = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(resultsFile)
    resultsFile = fullfile(here,'results','nested_cv_results_matlab.csv');
end
if nargin < 2
    referenceFile = fullfile(here,'results','reference_python','cv_results_all.csv');
end
R = readtable(resultsFile, 'TextType','string');
isDelta = startsWith(R.model,"DELTA");
A = R(~isDelta,:); D = R(isDelta,:);
outs = unique(A.outcome,'stable'); mets = unique(A.metric,'stable');

fprintf('\n=== Cross-validated performance, range across metrics ===\n');
fprintf('%-10s %-17s %-17s %-17s\n','outcome','tract16','tract4_tail','classical4');
for i = 1:numel(outs)
    fprintf('%-10s', outs(i));
    for m = ["tract16","tract4_tail","classical4"]
        v = A.cv(A.outcome == outs(i) & A.model == m);
        fprintf(' %5.2f to %5.2f   ', min(v), max(v));
    end
    fprintf('\n');
end

fprintf('\n=== Direction (cv point estimates) ===\n');
for m = ["tract16","tract4_tail"]
    wins = 0; tot = 0;
    for i = 1:numel(outs)
        for j = 1:numel(mets)
            s = A.outcome == outs(i) & A.metric == mets(j);
            t = A.cv(s & A.model == m); c = A.cv(s & A.model == "classical4");
            if ~isempty(t) && ~isempty(c), tot = tot + 1; wins = wins + (t > c); end
        end
    end
    fprintf('%s > classical4 in %d of %d metric-outcome combinations\n', m, wins, tot);
end

fprintf('\n=== Reliable difference (bootstrap 95%% interval above zero) ===\n');
for m = ["DELTA tract16 - classical4","DELTA tract4_tail - classical4"]
    d = D(D.model == m,:);
    fprintf('%s: %d of %d\n', m, sum(d.ci_lo > 0), height(d));
    for i = 1:numel(outs)
        k = d(d.outcome == outs(i) & d.ci_lo > 0,:);
        fprintf('   %-10s %s\n', outs(i), strjoin(k.metric', ', '));
    end
end

fprintf('\n=== Mean optimism (apparent minus cv) ===\n');
for ms = ["AUC","R2"]
    for m = ["tract16","tract4_tail","classical4"]
        s = A.measure == ms & A.model == m;
        fprintf('%-3s %-12s apparent %.3f  cv %.3f  optimism %.3f\n', ms, m, ...
            mean(A.apparent(s)), mean(A.cv(s)), mean(A.optimism(s)));
    end
end

if ~isempty(referenceFile) && isfile(referenceFile)
    P = readtable(referenceFile, 'TextType','string');
    J = innerjoin(R, P, 'Keys', {'metric','outcome','model'});
    fprintf('\n=== MATLAB vs Python reference (%d matched rows) ===\n', height(J));
    nMis = sum(J.N_R ~= J.N_P);
    fprintf('Rows with different N: %d (must be 0)\n', nMis);
    dcv = abs(J.cv_R - J.cv_P); dap = abs(J.apparent_R - J.apparent_P);
    fprintf('cv:       median |diff| %.3f, max %.3f\n', median(dcv,'omitnan'), max(dcv));
    fprintf('apparent: median |diff| %.3f, max %.3f\n', median(dap,'omitnan'), max(dap));
    [~, order] = sort(dcv,'descend');
    fprintf('Largest cv differences:\n');
    for k = order(1:min(8,numel(order)))'
        fprintf('   %-4s %-10s %-32s matlab %.3f  python %.3f\n', J.metric(k), J.outcome(k), ...
            J.model(k), J.cv_R(k), J.cv_P(k));
    end
end
end
