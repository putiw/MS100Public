function rev_p_correct()
% REV_P_CORRECT  Bonferroni correction with four independent families of 40 tests.
%
% Each file × sheet combination is treated as an independent family of 40 tests:
%   rev_lesionload_results.xlsx  "Results"  → 40 tests (no-demo)
%   rev_lesionload_results.xlsx  "Demo"     → 40 tests (+demo)
%   rev_qMRI_results.xlsx        "Results"  → 40 tests (no-demo)
%   rev_qMRI_results.xlsx        "Demo"     → 40 tests (+demo)
%
% Method: simple Bonferroni — p_adj = p * m, capped at 1.
% Consistent with the Bonferroni correction applied to the univariate
% correlation analysis (result252.m).
%
% Reads:
%   review/rev_lesionload_results.xlsx
%   review/rev_qMRI_results.xlsx
%
% Inserts a "p_adj" column immediately after "p" in each sheet and overwrites
% the files in place. Then calls rev_p_correct_bold.py to bold rows where
% p_adj < 0.05.
%
% Usage:
%   rev_p_correct

    clc;

    script_dir = fileparts(mfilename('fullpath'));

    files  = { ...
        fullfile(script_dir, 'rev_lesionload_results.xlsx'), ...
        fullfile(script_dir, 'rev_qMRI_results.xlsx') };
    sheets = {'Results', 'Demo'};

    for f = 1:numel(files)
        assert(isfile(files{f}), 'Missing: %s', files{f});
    end

    %% ── Load all four sheets ────────────────────────────────────────────
    tables    = cell(numel(files), numel(sheets));
    p_col_idx = zeros(numel(files), numel(sheets));

    for f = 1:numel(files)
        info = sheetnames(files{f});
        for s = 1:numel(sheets)
            assert(ismember(sheets{s}, info), ...
                'Sheet "%s" not found in %s', sheets{s}, files{f});
            T = readtable(files{f}, 'Sheet', sheets{s}, 'TextType', 'string');
            vn  = T.Properties.VariableNames;
            idx = find(strcmpi(vn, 'p'), 1);
            assert(~isempty(idx), 'No "p" column in %s [%s]', files{f}, sheets{s});
            assert(isnumeric(T.(vn{idx})), ...
                'p column not numeric in %s [%s]. Re-run rev_lesionload / rev_qMRI.', ...
                files{f}, sheets{s});
            tables{f, s}    = T;
            p_col_idx(f, s) = idx;
        end
    end

    %% ── Correct each sheet independently (40 tests each) ────────────────
    p_adj_store = cell(numel(files), numel(sheets));

    for f = 1:numel(files)
        for s = 1:numel(sheets)
            T    = tables{f, s};
            pidx = p_col_idx(f, s);
            np   = double(T.(T.Properties.VariableNames{pidx}));

            valid_mask = ~isnan(np);
            p_valid    = np(valid_mask);
            m          = sum(valid_mask);   % should be 40

            p_adj_all = nan(size(np));
            p_adj_all(valid_mask) = min(p_valid * m, 1);

            n_sig = sum(p_adj_all < 0.05, 'omitnan');
            [~, fname] = fileparts(files{f});
            fprintf('%-35s  [%s]  m=%d  sig=%d\n', fname, sheets{s}, m, n_sig);

            p_adj_store{f, s} = p_adj_all;
        end
    end

    %% ── Insert p_adj and write back ─────────────────────────────────────
    for f = 1:numel(files)
        for s = 1:numel(sheets)
            T    = tables{f, s};
            pidx = p_col_idx(f, s);
            p_adj_num = p_adj_store{f, s};

            % Drop any existing p_adj column to avoid duplication
            vn_right  = T.Properties.VariableNames((pidx+1):end);
            keep_right = ~strcmpi(vn_right, 'p_adj');

            left    = T(:, 1:pidx);
            right   = T(:, (pidx+1):end);
            right   = right(:, keep_right);
            p_adj_T = table(p_adj_num, 'VariableNames', {'p_adj'});

            T_new = [left, p_adj_T, right];
            writetable(T_new, files{f}, 'Sheet', sheets{s});
            fprintf('Written: %s [%s]  (%d rows)\n', files{f}, sheets{s}, height(T_new));
        end
    end

    %% ── Bold formatting ─────────────────────────────────────────────────
    py_script = fullfile(script_dir, 'rev_p_correct_bold.py');
    if ~isfile(py_script)
        warning('rev_p_correct_bold.py not found — skipping bold formatting.');
        return;
    end

    cmd = sprintf('python3 "%s" "%s" "%s"', py_script, files{1}, files{2});
    fprintf('\nRunning bold formatting:\n  %s\n', cmd);
    [status, out] = system(cmd);
    if status ~= 0
        warning('Python bold script failed (status %d):\n%s', status, out);
    else
        fprintf('%s\n', out);
    end
end
