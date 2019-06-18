curpath = fileparts(mfilename('fullpath'));
matfilename = fullfile(curpath, strcat('../../calculations/', clock2d.OUTPUT_MAT_FILENAME));
MEANVOL_FIGNAME = 'clock2d-meanvol.png';
SPECVOLS_FIGNAME = 'clock2d-specvols.png';
MEAN_VOL_BASE = 0.3;
MEAN_VOL_IMAGE_SCALING = 4;
SPECTRAL_VOL_IMAGE_SCALING = 3;

fprintf('Loading %s.max(mat\n', matfilename);
load(matfilename, 'specvols')
mean_vol = ((specvols(:,:,1)/sqrt(clock2d.n))-MEAN_VOL_BASE)*MEAN_VOL_IMAGE_SCALING;
specvols_without_mean = specvols(:,:,2:end)*SPECTRAL_VOL_IMAGE_SCALING/sqrt(clock2d.n) + 0.5;

figpath = fullfile(curpath, '../../figures/');

tic;

meanvol_filename = fullfile(figpath, MEANVOL_FIGNAME); 
fprintf('Saving to %s\n', meanvol_filename);
imwrite(mean_vol, meanvol_filename);

specvol_filename = fullfile(figpath, SPECVOLS_FIGNAME);
fprintf('Saving to %s\n', specvol_filename);
assert(mod(clock2d.r,2)==1);
imwrite(imtile(specvols_without_mean(:,:,1:clock2d.NUM_SPECVOLS_TO_DISPLAY-1), 'GridSize', [(clock2d.NUM_SPECVOLS_TO_DISPLAY-1)/2 2])', specvol_filename);

toc;
%set(gcf,'Units','centimeters');
%screenposition = get(gcf,'Position');
%set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], ...
%    'PaperSize', [screenposition(3:4)]);
%saveas(gcf, strcat(figfilename,'.pdf'));
%imwrite(tiled, strcat(figfilename,'.png'));
%print('-deps', figfilename);
