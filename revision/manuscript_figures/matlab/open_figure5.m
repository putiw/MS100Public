function fig = open_figure5()
% Open an editable MATLAB figure from frozen aggregates without exporting.
revisionDir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(revisionDir,'code'));
fig=figure5_publication([],[],false);
end
