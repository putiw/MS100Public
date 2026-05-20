function table1_verify()
% TABLE1_VERIFY  Compute Table 1 summary statistics from clinicalScore.xlsx
%
% Reports per group (Controls, RR, PP, SP, All MS, Overall) and overall:
%   N (% Female)
%   Age (mean ± SD)
%   Disease duration (mean ± SD)
%   EDSS (median [IQR])
%   On DMT: N (%)
%   Therapy 1st line: N (%)
%   Therapy 2nd line: N (%)
%   T25FW (mean ± SD)
%   9HPT-D (mean ± SD)
%   9HPT-ND (mean ± SD)
%   SDMT Total (mean ± SD)
%   MSPro (median [IQR])
%
% Group coding in data: 0=Controls, 1=RR, 2=PP, 3=SP

cfg = load_config();
data_file = fullfile(cfg.bidsDir, cfg.statsDir, 'clinicalScore.xlsx');

fprintf('Loading data from:\n  %s\n\n', data_file);
T = readtable(data_file);

% ── Column mapping (verify against header) ────────────────────────────────
fprintf('Columns found: %s\n\n', strjoin(T.Properties.VariableNames, ', '));

% Group definitions
group_codes = {0, 1, 2, 3};
group_names = {'Controls', 'RR', 'PP', 'SP'};

% Build index for each group and useful composites
idx = struct();
for gi = 1:4
    idx.(group_names{gi}) = T.Group == group_codes{gi};
end
idx.AllMS    = T.Group == 1 | T.Group == 2 | T.Group == 3;
idx.Overall  = true(height(T), 1);

col_order = {'Controls', 'RR', 'PP', 'SP', 'AllMS', 'Overall'};

% ── Helper functions ───────────────────────────────────────────────────────
mean_sd = @(v) sprintf('%.1f ± %.1f', mean(v,'omitnan'), std(v,'omitnan'));
med_iqr = @(v) sprintf('%.1f [%.1f–%.1f]', ...
    median(v,'omitnan'), quantile(v(isfinite(v)), 0.25), quantile(v(isfinite(v)), 0.75));
n_pct   = @(mask_num, mask_denom) sprintf('%d (%.1f%%)', ...
    sum(mask_num,'omitnan'), 100*sum(mask_num,'omitnan')/max(sum(mask_denom,'omitnan'),1));

% ── Print table ────────────────────────────────────────────────────────────
hdr = sprintf('%-30s', 'Variable');
for c = col_order
    hdr = [hdr sprintf('%20s', c{1})];
end
sep = repmat('-', 1, length(hdr));
fprintf('%s\n%s\n%s\n', sep, hdr, sep);

% ── N ─────────────────────────────────────────────────────────────────────
row = sprintf('%-30s', 'N');
for c = col_order
    n = sum(idx.(c{1}));
    row = [row sprintf('%20d', n)];
end
fprintf('%s\n', row);

% ── % Female ──────────────────────────────────────────────────────────────
row = sprintf('%-30s', '% Female');
for c = col_order
    mask = idx.(c{1});
    is_female = mask & strcmp(T.Gender, 'F');
    row = [row sprintf('%20s', n_pct(is_female, mask))];
end
fprintf('%s\n', row);

% ── Age ───────────────────────────────────────────────────────────────────
row = sprintf('%-30s', 'Age (mean ± SD)');
for c = col_order
    v = T.Age(idx.(c{1}));
    row = [row sprintf('%20s', mean_sd(v))];
end
fprintf('%s\n', row);

% ── Disease duration ──────────────────────────────────────────────────────
row = sprintf('%-30s', 'Disease duration (mean ± SD)');
for c = col_order
    v = T.DurationOfDisease(idx.(c{1}));
    if strcmp(c{1}, 'Controls')
        row = [row sprintf('%20s', 'N/A')];
    else
        row = [row sprintf('%20s', mean_sd(v))];
    end
end
fprintf('%s\n', row);

% ── EDSS ─────────────────────────────────────────────────────────────────
row = sprintf('%-30s', 'EDSS (median [IQR])');
for c = col_order
    v = T.EDSS(idx.(c{1}));
    v = v(isfinite(v));
    if strcmp(c{1}, 'Controls')
        row = [row sprintf('%20s', '0 [0–0]')];
    else
        row = [row sprintf('%20s', med_iqr(v))];
    end
end
fprintf('%s\n', row);

% ── On DMT ───────────────────────────────────────────────────────────────
row = sprintf('%-30s', 'On DMT: N (%)');
for c = col_order
    mask = idx.(c{1});
    if strcmp(c{1}, 'Controls')
        row = [row sprintf('%20s', 'N/A')];
    else
        on_dmt = mask & T.AdministrationOfTherapy == 1;
        row = [row sprintf('%20s', n_pct(on_dmt, mask))];
    end
end
fprintf('%s\n', row);

% ── Therapy 1st line ─────────────────────────────────────────────────────
row = sprintf('%-30s', 'Therapy 1st line: N (%)');
for c = col_order
    mask = idx.(c{1});
    if strcmp(c{1}, 'Controls')
        row = [row sprintf('%20s', 'N/A')];
    else
        line1 = mask & ~isnan(T.DrugLine) & T.DrugLine == 1;
        row = [row sprintf('%20s', n_pct(line1, mask))];
    end
end
fprintf('%s\n', row);

% ── Therapy 2nd line ─────────────────────────────────────────────────────
row = sprintf('%-30s', 'Therapy 2nd line: N (%)');
for c = col_order
    mask = idx.(c{1});
    if strcmp(c{1}, 'Controls')
        row = [row sprintf('%20s', 'N/A')];
    else
        line2 = mask & ~isnan(T.DrugLine) & T.DrugLine == 2;
        row = [row sprintf('%20s', n_pct(line2, mask))];
    end
end
fprintf('%s\n', row);

% ── T25FW ────────────────────────────────────────────────────────────────
row = sprintf('%-30s', 'T25FW (mean ± SD)');
for c = col_order
    v = T.T25FW(idx.(c{1}));
    v = v(isfinite(v));
    row = [row sprintf('%20s', mean_sd(v))];
end
fprintf('%s\n', row);

% ── 9HPT-D ───────────────────────────────────────────────────────────────
row = sprintf('%-30s', '9HPT-D (mean ± SD)');
for c = col_order
    v = T.x9HPTD(idx.(c{1}));
    v = v(isfinite(v));
    row = [row sprintf('%20s', mean_sd(v))];
end
fprintf('%s\n', row);

% ── 9HPT-ND ──────────────────────────────────────────────────────────────
row = sprintf('%-30s', '9HPT-ND (mean ± SD)');
for c = col_order
    v = T.x9HPTND(idx.(c{1}));
    v = v(isfinite(v));
    row = [row sprintf('%20s', mean_sd(v))];
end
fprintf('%s\n', row);

% ── SDMT Total ───────────────────────────────────────────────────────────
row = sprintf('%-30s', 'SDMT Total (mean ± SD)');
for c = col_order
    v = T.SDMTtotal(idx.(c{1}));
    v = v(isfinite(v));
    row = [row sprintf('%20s', mean_sd(v))];
end
fprintf('%s\n', row);

% ── MSPro ────────────────────────────────────────────────────────────────
% MSPro is ordinal (1/2/3): report median [IQR] for MS groups, N per category
row = sprintf('%-30s', 'MSPro (median [IQR])');
for c = col_order
    v = T.MSPro(idx.(c{1}));
    v = v(isfinite(v));
    if strcmp(c{1}, 'Controls') || isempty(v)
        row = [row sprintf('%20s', 'N/A')];
    else
        row = [row sprintf('%20s', med_iqr(v))];
    end
end
fprintf('%s\n', row);

% MSPro breakdown by category
for cat = 1:3
    row = sprintf('%-30s', sprintf('  MSPro = %d: N (%%)', cat));
    for c = col_order
        mask = idx.(c{1});
        if strcmp(c{1}, 'Controls')
            row = [row sprintf('%20s', 'N/A')];
        else
            v = T.MSPro(mask);
            has_mspro = mask & ~isnan(T.MSPro);
            cat_mask  = mask & ~isnan(T.MSPro) & T.MSPro == cat;
            row = [row sprintf('%20s', n_pct(cat_mask, has_mspro))];
        end
    end
    fprintf('%s\n', row);
end

fprintf('%s\n\n', sep);

% ── Missing data summary ──────────────────────────────────────────────────
fprintf('=== MISSING DATA SUMMARY (MS patients only) ===\n');
ms_idx = idx.AllMS;
vars_check = {'T25FW','x9HPTD','x9HPTND','SDMTtotal','MSPro','EDSS','DurationOfDisease'};
for vi = 1:numel(vars_check)
    v = T.(vars_check{vi})(ms_idx);
    n_miss = sum(isnan(v));
    fprintf('  %-25s: %d missing / %d total\n', vars_check{vi}, n_miss, sum(ms_idx));
end

end
