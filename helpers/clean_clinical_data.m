function data = clean_clinical_data(data, clinical_scores)
    % Clean and preprocess clinical data
    %
    % Input:
    %   data: Table containing the data
    %   clinical_scores: Cell array of clinical score names
    %
    % Output:
    %   data: Table with cleaned clinical data
    
    % log-transform T25FW first
    if ismember('T25FW', data.Properties.VariableNames)
        data.T25FW = log(data.T25FW);
    end
    
    % winsorize each score
    for ii = 1:numel(clinical_scores)
        fld = clinical_scores{ii};
        if ismember(fld, data.Properties.VariableNames)
            col = data.(fld);
            m   = mean(col,'omitnan');
            sd  = std(col,'omitnan');
            hi  = m + 2.5*sd;
            lo  = m - 2.5*sd;
            % clamp extreme values
            data.(fld)(col > hi) = hi;
            data.(fld)(col < lo) = lo;
        end
    end
end 