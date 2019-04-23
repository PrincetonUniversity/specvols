MAT_FILENAME = 'clock2d_diffvols';
MEANVOL_FIGNAME = 'clock2d_meanvol.png';
DIFFVOLS_FIGNAME = 'clock2d_diffvols.png';

curpath = fileparts(mfilename('fullpath'));
matfilename = fullfile(curpath, strcat('../../calculations/', MAT_FILENAME));
fprintf('Loading %s.mat\n', matfilename);
load(matfilename, 'diffvols')
mean_vol = diffvols(:,:,1);
SCALE = 32;
diffvols_without_meanvol = diffvols(:,:,2:end)/SCALE + 0.5;

figpath = fullfile(curpath, '../../figures/');

tic;

meanvol_filename = fullfile(figpath, MEANVOL_FIGNAME); 
fprintf('Saving to %s\n', meanvol_filename);
imwrite(mean_vol, meanvol_filename);

diffvol_filename = fullfile(figpath, DIFFVOLS_FIGNAME);
fprintf('Saving to %s\n', diffvol_filename);
imwrite(imtile(diffvols_without_meanvol), diffvol_filename);

toc;
%set(gcf,'Units','centimeters');
%screenposition = get(gcf,'Position');
%set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], ...
%    'PaperSize', [screenposition(3:4)]);
%saveas(gcf, strcat(figfilename,'.pdf'));
%imwrite(tiled, strcat(figfilename,'.png'));
%print('-deps', figfilename);
