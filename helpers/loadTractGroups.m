function groupDefs = loadTractGroups()
% LOADTRACTGROUPS Loads all tract group definitions from configs/tract_groups directory
%
% Returns:
%   groupDefs - Structure where each field is a group name and contains a cell
%               array of tract names belonging to that group

    scriptPath = fileparts(fileparts(mfilename('fullpath'))); % Go up one level to Clean/
    groupsDir = fullfile(scriptPath, 'configs', 'tract_groups');
    
    % List all tract group files
    groupFiles = dir(fullfile(groupsDir, '*.txt'));
    groupDefs = struct();
    
    for i = 1:length(groupFiles)
        filePath = fullfile(groupsDir, groupFiles(i).name);
        groupDefs = parseTractGroupFile(filePath, groupDefs);
    end
end

