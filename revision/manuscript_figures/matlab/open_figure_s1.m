function fig = open_figure_s1()
% Open an editable MATLAB figure from frozen aggregates without exporting.
revisionDir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(revisionDir,'code'));
[~,fig]=supplementary_figure_s1_correlations_sdmt_msfc('doExport',false,'keepOpen',true);
end
