function fig = open_manuscript_artwork(number,pdftoppm)
% Display artwork as an editable MATLAB image object using a disposable render.
name=sprintf('Figure%d.pdf',number);
package=fileparts(fileparts(fileparts(mfilename('fullpath'))));
source=fullfile(package,'formatting','originals',name);
if ~isfile(source)
    workspace=fileparts(fileparts(fileparts(package)));
    source=fullfile(workspace,'final_figures_20261106',name);
end
assert(isfile(source),'Missing local artwork for Figure %d.',number);
scratch=tempname;
mkdir(scratch);
cleanup=onCleanup(@()rmdir(scratch,'s')); %#ok<NASGU>
preview=fullfile(scratch,'artwork');
command=strjoin({shell_quote(pdftoppm),'-r 200 -png -singlefile',shell_quote(source),shell_quote(preview)},' ');
[status,output]=system(command);
assert(status==0,'Cannot open Figure %d. Check that pdftoppm is available: %s',number,output);
pixels=imread([preview '.png']);
fig=figure('Name',sprintf('Figure %d',number));
ax=axes(fig);
image(ax,pixels);
axis(ax,'image');axis(ax,'off');
fig.UserData=source;
end

function quoted=shell_quote(value)
apostrophe=char(39);
quoted=[apostrophe strrep(char(value),apostrophe,[apostrophe '"' apostrophe '"' apostrophe]) apostrophe];
end
