function data = load_and_merge_data(cfg, varargin)
    % Load and merge data from various sources
    %
    % Input:
    %   cfg: Configuration structure
    %   varargin: Optional parameters
    %       'data_types': Cell array of data types to load ('clinical', 'tract', 'ico', 'T1', 'MD', 'MTR', 'FA', 'AD', 'RD')
    %       'drop_controls': Boolean, whether to drop control subjects (default: false)
    %
    % Output:
    %   data: Merged table containing all requested data
    
    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'data_types', {'clinical', 'tract'}, @iscell);
    addParameter(p, 'drop_controls', false, @islogical);
    parse(p, varargin{:});
    
    baseDir = fullfile(cfg.bidsDir, cfg.statsDir);
    
    % Initialize with clinical data
    if ismember('clinical', p.Results.data_types)
        data = readtable(fullfile(baseDir, 'clinicalScore.xlsx'));
    else
        error('Clinical data must be included');
    end
    
    % Load tract-based lesion load data
    if ismember('tract', p.Results.data_types)
        Ttract = readtable(fullfile(baseDir, 'TractLesionLoad.xlsx'));
        data = outerjoin(data, Ttract, 'Keys', 'SubjectID', 'MergeKeys', true);
    end
    
    % Load Icometrix data
    if ismember('ico', p.Results.data_types)
        Tico = readtable(fullfile(baseDir, 'icometrixLesionStats.xlsx'));
        data = outerjoin(data, Tico, 'Keys', 'SubjectID', 'MergeKeys', true);
    end
    
    % Load MRI metrics
    if ismember('T1', p.Results.data_types)
        T1 = readtable(fullfile(baseDir, 'GroupTractT1.xlsx'));
        data = outerjoin(data, T1, 'Keys', 'SubjectID', 'MergeKeys', true);
    end
    if ismember('MD', p.Results.data_types)
        MD = readtable(fullfile(baseDir, 'GroupTractMD.xlsx'));
        data = outerjoin(data, MD, 'Keys', 'SubjectID', 'MergeKeys', true);
    end
    if ismember('MTR', p.Results.data_types)
        MTR = readtable(fullfile(baseDir, 'GroupTractMTR.xlsx'));
        data = outerjoin(data, MTR, 'Keys', 'SubjectID', 'MergeKeys', true);
    end
    if ismember('FA', p.Results.data_types)
        FA = readtable(fullfile(baseDir, 'GroupTractFA.xlsx'));
        data = outerjoin(data, FA, 'Keys', 'SubjectID', 'MergeKeys', true);
    end
    if ismember('AD', p.Results.data_types)
        AD = readtable(fullfile(baseDir, 'GroupTractAD.xlsx'));
        data = outerjoin(data, AD, 'Keys', 'SubjectID', 'MergeKeys', true);
    end
    if ismember('RD', p.Results.data_types)
        RD = readtable(fullfile(baseDir, 'GroupTractRD.xlsx'));
        data = outerjoin(data, RD, 'Keys', 'SubjectID', 'MergeKeys', true);
    end
    
    % Drop control subjects if requested
    if p.Results.drop_controls
        data(startsWith(data.SubjectID, 'sub-C'), :) = [];
    end
end 