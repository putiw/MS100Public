function groupDefs = parseTractGroupFile(filePath, groupDefs)
% PARSETRACTGROUPFILE Parses a single tract group file
%
% Usage:
%   groupDefs = parseTractGroupFile(filePath)
%   groupDefs = parseTractGroupFile(filePath, existingGroupDefs)
%
% Input:
%   filePath - Path to the tract group file
%   groupDefs - (Optional) Existing group definitions structure to append to
%
% Returns:
%   groupDefs - Updated structure containing tract group definitions
%
% File Format:
%   [GroupName]
%   TractName1
%   TractName2
%   ...

    % Initialize empty struct if not provided
    if nargin < 2
        groupDefs = struct();
    end
    
    fid = fopen(filePath, 'r');
    if fid == -1
        error('Could not open tract group file: %s', filePath);
    end
    
    currentGroup = '';
    
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        
        % Skip empty lines and comments
        if isempty(line) || startsWith(line, '#')
            continue;
        end
        
        % Check if this is a group header
        if startsWith(line, '[') && endsWith(line, ']')
            currentGroup = line(2:end-1);
            groupDefs.(currentGroup) = {};
            continue;
        end
        
        % Add tract to current group
        if ~isempty(currentGroup)
            groupDefs.(currentGroup){end+1} = line;
        end
    end
    
    fclose(fid);
end