function config = load_config()
%  Loads and validates configuration from configs/config.json
%
% Returns:
%   config - Structure containing configuration settings

    % Get path to config file
    scriptPath = fileparts(mfilename('fullpath'));
    configFile = fullfile(scriptPath, '..', 'configs', 'config.json');
    
    % Read and parse JSON config
    fid = fopen(configFile, 'r');
    if fid == -1
        error('Could not open config file: %s', configFile);
    end
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    
    config = jsondecode(str);
    
    % Validate required fields
    requiredFields = {'bidsDir', 'derivDir', 'tractTypes'};
    for i = 1:length(requiredFields)
        if ~isfield(config, requiredFields{i})
            error('Missing required field in config.json: %s', requiredFields{i});
        end
    end
    
    % Add default values for optional fields
    if ~isfield(config, 'excelOutputs')
        fprintf('Warning: excelOutputs not found in config.json, adding defaults\n');
        config.excelOutputs = struct(...
            'individual', 'lesion_load_individual.xlsx', ...
            'group', 'lesion_load_group.xlsx', ...
            'templateGroup', 'lesion_load_template_group.xlsx' ...
        );
    end
end 