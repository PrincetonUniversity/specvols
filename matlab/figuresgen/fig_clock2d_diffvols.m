MAT_FILENAME = 'clock2d_diffvols';
MEANVOL_FIGNAME = 'clock2d_meanvol';
DIFFVOLS_FIGNAME = 'clock2d_diffvols';

matfilename = fullfile(fileparts(mfilename('fullpath')), strcat('../../calculations/', MAT_FILENAME));
fprintf('Loading %s.mat\n', matfilename);
load(matfilename, 'diffvols')
diffvols_without_meanvol = diffvols(:,:,2:end);
SCALE = 20;
montage(diffvols(:,:,2:end)/SCALE+0.5);

figfilename = fullfile(fileparts(mfilename('fullpath')), strcat('../../figures/', DIFFVOLS_FIGNAME));
fprintf('Saving to %s\n', figfilename);
print(figfilename,'-dpng');