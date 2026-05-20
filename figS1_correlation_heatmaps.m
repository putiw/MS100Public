function figS1_correlation_heatmaps()
% FIGS1_CORRELATION_HEATMAPS  Supplementary Figure S1: partial correlations
% for whole-brain lesion load metrics.
%
% Four rev_corr calls form one family of 130 tests:
%   WB (2 rows × 5 outcomes = 10)
%   Tract-based LN / LV / Lnorm (4 groups × 5 outcomes × 3 = 60)
%   Classical LN / LV / Lnorm   (4 regions × 5 outcomes × 3 = 60)
%
% Correction: Bonferroni applied jointly over all 130 tests.
% Figures are drawn with the Bonferroni-corrected p-values (* = p_adj < 0.05).

close all;

%% ── Collect raw data without plotting ────────────────────────────────────
r_wb  = rev_corr('WB',     'plot', false, 'Whole Brain Lesion Load');
r_num = rev_corr('number', 'groupOnly', true, 'plot', false, 'Lesion Number');
r_vol = rev_corr('volume', 'groupOnly', true, 'plot', false, 'Lesion Volume');
r_nor = rev_corr('norm',   'groupOnly', true, 'plot', false, 'Normalized Lesion Volume');

all_results = [r_wb(:); r_num(:); r_vol(:); r_nor(:)];   % 7 heatmaps total

%% ── Pool all raw p-values ────────────────────────────────────────────────
all_p  = [];
sizes  = zeros(numel(all_results), 1);
for i  = 1:numel(all_results)
    pv       = all_results(i).p_raw(:);
    sizes(i) = numel(pv);
    all_p    = [all_p; pv]; %#ok<AGROW>
end
offsets = [0; cumsum(sizes(1:end-1))];

fprintf('Total tests: %d  (non-NaN: %d)\n', numel(all_p), sum(~isnan(all_p)));

%% ── Bonferroni correction over all 130 tests ─────────────────────────────
p_bonf_all = bonferroni(all_p);
fprintf('Significant after Bonferroni (p_adj < 0.05): %d / %d\n', ...
    sum(p_bonf_all < 0.05, 'omitnan'), sum(~isnan(all_p)));

%% ── Print Whole Brain results ────────────────────────────────────────────
wb         = all_results(1);   % r_wb always first
idx_wb     = offsets(1) + (1:sizes(1));
p_wb       = reshape(p_bonf_all(idx_wb), size(wb.p_raw));

fprintf('\n%s\n  Whole Brain Lesion Load — r and Bonferroni-adjusted p\n%s\n', ...
    repmat('─',1,65), repmat('─',1,65));
fprintf('  %-22s', '');
for j = 1:numel(wb.clinical_labels)
    fprintf('  %-14s', wb.clinical_labels{j});
end
fprintf('\n');
for k = 1:numel(wb.display_names)
    fprintf('  %-22s', wb.display_names{k});
    for j = 1:numel(wb.clinical_labels)
        r    = wb.corr_matrix(k,j);
        padj = p_wb(k,j);
        if isnan(r)
            fprintf('  %-14s', 'n/a');
        else
            sig = '';
            if padj < 0.05, sig = '*'; end
            if padj < 0.001
                fprintf('  r=%.3f p<0.001%-2s', r, sig);
            else
                fprintf('  r=%.3f p=%.3f%-2s', r, padj, sig);
            end
        end
    end
    fprintf('\n');
end
fprintf('%s\n  * significant after Bonferroni correction (p_adj < 0.05)\n%s\n\n', ...
    repmat('─',1,65), repmat('─',1,65));

%% ── Replot with Bonferroni-corrected p-values ────────────────────────────
for i = 1:numel(all_results)
    idx         = offsets(i) + (1:sizes(i));
    p_corrected = reshape(p_bonf_all(idx), size(all_results(i).p_raw));
    rev_corr_replot(all_results(i), p_corrected);
end


%% ── Bonferroni helper ────────────────────────────────────────────────────
function p_adj = bonferroni(p_values)
% Bonferroni correction: p_adj = p * n_valid, capped at 1. NaNs are preserved.
    p_flat       = p_values(:);
    n            = sum(~isnan(p_flat));
    p_adj        = min(p_flat * n, 1);
    p_adj(isnan(p_flat)) = NaN;
    p_adj        = reshape(p_adj, size(p_values));
end
