function rev_corr_replot(result, p_corrected)
% REV_CORR_REPLOT  Replot a rev_corr heatmap with externally-corrected p-values.
%
% Inputs:
%   result      — one element from the struct array returned by rev_corr()
%   p_corrected — matrix of corrected p-values, same size as result.p_raw

    corr_matrix     = result.corr_matrix;
    display_names   = result.display_names;
    clinical_labels = result.clinical_labels;
    title_str       = result.title;

    num_tracts = size(corr_matrix, 1);
    num_scores = size(corr_matrix, 2);

    fig_w_in    = 3.5;
    cell_h_in   = 0.20;
    margin_h_in = 1.1;
    fig_h_in    = margin_h_in + num_tracts * cell_h_in;

    figure('Units','inches', 'Position',[1 1 fig_w_in fig_h_in], ...
           'PaperUnits','inches', 'PaperSize',[fig_w_in fig_h_in], ...
           'DefaultAxesFontName','Arial', 'DefaultTextFontName','Arial');

    imagesc(corr_matrix);
    colorbar;
    colormap(flip(hot));
    caxis([0 1]);

    title([title_str ' Correlations'], 'FontSize',7, 'FontWeight','bold');

    ax_fs = 6;
    set(gca, 'XTick',1:num_scores, 'XTickLabel',clinical_labels, 'FontSize',ax_fs);
    set(gca, 'YTick',1:num_tracts, 'YTickLabel',display_names,   'FontSize',ax_fs);
    xtickangle(30);

    cell_fs = 6;
    for k = 1:num_tracts
        for j = 1:num_scores
            if ~isnan(corr_matrix(k,j))
                sig = p_corrected(k,j) < 0.05;
                txt = sprintf('%.2f', corr_matrix(k,j));
                if sig, txt = [txt '*']; end %#ok<AGROW>
                fw = 'normal'; if sig, fw = 'bold'; end
                fc = 'k';      if corr_matrix(k,j) > 0.6, fc = 'w'; end
                text(j, k, txt, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'FontSize',cell_fs, 'FontWeight',fw, 'Color',fc);
            end
        end
    end
end
