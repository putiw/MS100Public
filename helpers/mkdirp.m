function mkdirp(dirPath)
% MKDIRP Creates directory recursively (like mkdir -p in Unix)
%
% Usage:
%   mkdirp('/path/to/new/directory')
%
% Args:
%   dirPath: Path to create
%
% Creates all directories in the path if they don't exist.
% Does nothing if the directory already exists.

    if ~exist(dirPath, 'dir')
        mkdir(dirPath);
    end
end 