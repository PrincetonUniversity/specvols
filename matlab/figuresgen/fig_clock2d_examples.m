curpath = fileparts(mfilename('fullpath'));
matfilename = fullfile(curpath, strcat('../../calculations/', clock2d.OUTPUT_MAT_FILENAME));
fprintf('Loading %s.mat\n', matfilename);
load(matfilename, 'specvols', 'clock_examples_clean', 'clock_examples_noisy', 'clock_examples_reconstructed');

assert(isequal(size(clock_examples_clean), size(clock_examples_noisy)));
assert(isequal(size(clock_examples_clean), size(clock_examples_reconstructed)));

curpath = fileparts(mfilename('fullpath'));
figpath = fullfile(curpath, '../../figures/');
for i=1:size(clock_examples_clean,1)
    fn = fullfile(figpath, strcat('clock2d-example-clean-', int2str(i), '.png'));
    fprintf('Saving to %s\n', fn);
    imwrite(squeeze(clock_examples_clean(i,:,:)), fn);

    fn = fullfile(figpath, strcat('clock2d-example-noisy-', int2str(i), '.png'));
    fprintf('Saving to %s\n', fn);
    imwrite(squeeze(clock_examples_noisy(i,:,:)), fn);

    fn = fullfile(figpath, strcat('clock2d-example-reconstructed-', int2str(i), '.png'));
    fprintf('Saving to %s\n', fn);
    imwrite((squeeze(clock_examples_reconstructed(i,:,:))), fn);
end
