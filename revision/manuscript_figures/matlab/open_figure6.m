function fig = open_figure6()
% Open an editable MATLAB figure from frozen aggregates without exporting.
revisionDir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(revisionDir,'code'));
[~,fig]=figure6_r2_comparison_sdmt_msfc_noddi('doExport',false,'keepOpen',true,'variants',1);
end
