function [corr_matrix, p_matrix] = calculate_correlations(data, metrics, clinical_scores, varargin)
    % Calculate correlations between metrics and clinical scores
    %
    % Input:
    %   data: Table containing the data
    %   metrics: Cell array of metric names
    %   clinical_scores: Cell array of clinical score names
    %   varargin: Optional parameters
    %       'skip_zeros': Boolean, whether to skip zero values (default: false)
    %       'skip_metric_zeros': Boolean, whether to skip zero values in metrics (default: false)
    %       'skip_clinical_zeros': Boolean, whether to skip zero values in clinical scores (default: false)
    %
    % Output:
    %   corr_matrix: Matrix of absolute correlation values
    %   p_matrix: Matrix of p-values
    
    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'skip_zeros', false, @islogical);
    addParameter(p, 'skip_metric_zeros', false, @islogical);
    addParameter(p, 'skip_clinical_zeros', false, @islogical);
    parse(p, varargin{:});
    
    % Initialize correlation matrix
    num_metrics = length(metrics);
    num_scores = length(clinical_scores);
    corr_matrix = zeros(num_metrics, num_scores);
    p_matrix = zeros(num_metrics, num_scores);
    
    % Calculate correlations for each metric-score pair
    for k = 1:num_metrics
        metric = metrics{k};
        for j = 1:num_scores
            clinical_score = clinical_scores{j};
            
            % Remove rows with missing values
            valid_rows = ~isnan(data.(clinical_score)) & ...
                        ~isnan(data.(metric));
            
            % Handle zero values based on parameters
            if p.Results.skip_zeros || p.Results.skip_metric_zeros
                valid_rows = valid_rows & data.(metric) ~= 0;
            end
            if p.Results.skip_zeros || p.Results.skip_clinical_zeros
                valid_rows = valid_rows & data.(clinical_score) ~= 0;
            end
            
            analysis_data = data(valid_rows, :);
            
            % Skip if not enough data points
            if height(analysis_data) < 3
                corr_matrix(k, j) = NaN;
                p_matrix(k, j) = NaN;
                continue;
            end
            
            % Calculate correlation
            [r, p_val] = corr(analysis_data.(metric), analysis_data.(clinical_score));
            corr_matrix(k, j) = abs(r);  % Use absolute value
            p_matrix(k, j) = p_val;
        end
    end
end 