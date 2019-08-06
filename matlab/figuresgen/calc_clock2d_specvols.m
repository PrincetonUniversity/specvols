
L = clock2d.L;
n = clock2d.n;
r = clock2d.r;
GAUSSIAN_BLUR_SIGMA = clock2d.GAUSSIAN_BLUR_SIGMA;
NOISE_SIGMA = clock2d.NOISE_SIGMA;
GAUSSIAN_KERNEL_SIGMA = 0.5*sqrt(L*L)*NOISE_SIGMA;
OUTPUT_MAT_FILENAME = clock2d.OUTPUT_MAT_FILENAME;
CLOCK_EXAMPLE_ANGLES = [6, 4.5, 3, 1];

assert(mod(L,2)==0);
H = L/2;
[X,Y] = meshgrid((-H:H-1)/H,(-H:H-1)/H);
clockface2d = X.^2 + Y.^2 <= (0.8^2+0.1^2);
clockstand_coords = [0 0; 0.5 -1; -0.5 -1; 0 0];
clockstand2d = poly2mask(clockstand_coords(:,1)*(-H)+H+1, clockstand_coords(:,2)*(-H)+H+1, L,L);
clock_2D_body = clockface2d|clockstand2d;

angles = 2*pi*rand(n,1);
angles = sort(angles);

tic;
fprintf('Generating clock face vectors... ');
clocks_clean = zeros(n, L^2);
clocks_noisy = zeros(n, L^2);
for i=1:n
    clock_2D_hand = gen_clock_hand_2d(L, angles(i), GAUSSIAN_BLUR_SIGMA);
    clock_2D = 0.4*clock_2D_body + clock_2D_hand;
    flattenedclock = reshape(clock_2D, 1, L^2);
    clocks_clean(i,:) = flattenedclock;
    clocks_noisy(i,:) = flattenedclock + normrnd(0, NOISE_SIGMA, size(flattenedclock));
end
fprintf('Finished.\n');
toc;

tic;
fprintf('Forming Laplacian... ');
W = graph_gaussian_kernel(clocks_noisy, GAUSSIAN_KERNEL_SIGMA);
lap = laplacian(W, 'normalized');
fprintf('Done.');
toc;

tic;
fprintf('Computing eigenvectors... ');
[laplacian_eigvecs,D] = eigs(lap, r, 'smallestabs');
fprintf('Done.');
toc;

% Eigenvectors have arbitrary sign. We want the first eigenvector
% (that is constant for the combinatorial laplacian)
% to be positive.
if sum(laplacian_eigvecs(:,1)) < 0
    laplacian_eigvecs(:,1) = -laplacian_eigvecs(:,1);
end

tic;
fprintf('Solving least-squares problem... ');
specvols = mldivide(laplacian_eigvecs, clocks_noisy);
fprintf('Done.');
toc;

specvols = reshape(specvols, r, L, L);
specvols = permute(specvols, [2 3 1]);

clock_examples_clean = zeros(length(CLOCK_EXAMPLE_ANGLES), L, L);
clock_examples_noisy = zeros(length(CLOCK_EXAMPLE_ANGLES), L, L);
clock_examples_reconstructed = zeros(length(CLOCK_EXAMPLE_ANGLES), L, L);

for i=1:length(CLOCK_EXAMPLE_ANGLES)
    angle = CLOCK_EXAMPLE_ANGLES(i);
    index = round(n*angle/(2*pi));
    clock_examples_clean(i,:,:) = reshape(clocks_clean(index,:), [L L]);
    clock_examples_noisy(i,:,:) = reshape(clocks_noisy(index,:), [L L]);
    reconstructed = reshape(specvols, [L*L r]) * laplacian_eigvecs(index, :)';
    clock_examples_reconstructed(i,:,:) = reshape(reconstructed, [L L]);
end


filename = fullfile(fileparts(mfilename('fullpath')), strcat('../../calculations/', OUTPUT_MAT_FILENAME));
fprintf('Saving to %s.mat\n', filename);
save(filename, 'L', 'n', 'r', 'GAUSSIAN_KERNEL_SIGMA', 'laplacian_eigvecs', 'specvols', 'clock_examples_clean', 'clock_examples_noisy','clock_examples_reconstructed');
