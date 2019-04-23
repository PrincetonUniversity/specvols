n = 2^15;
L = 127;
r = 17;
GAUSSIAN_KERNEL_SIGMA = 200;
OUTPUT_MAT_FILENAME = 'clock2d_diffvols';

angles = 2*pi*rand(n,1);
angles = sort(angles);

clocks = zeros(n,L^2);
for i=1:n
    clock_2D = gen_clock_hand_2d(L, angles(i));
    flattenedclock = reshape(clock_2D, 1, L^2);
    clocks(i,:) = flattenedclock;
end

W = graph_gaussian_kernel(clocks, GAUSSIAN_KERNEL_SIGMA);

lap = laplacian(W, 'combinatorial');

[laplacian_eigvecs,D] = eigs(lap, r, 'smallestabs');

diffvols = mldivide(laplacian_eigvecs, clocks);
diffvols = reshape(diffvols,r,L,L);
diffvols = permute(diffvols, [2 3 1]);

filename = fullfile(fileparts(mfilename('fullpath')), strcat('../../calculations/', OUTPUT_MAT_FILENAME));
fprintf('Saving to %s.mat\n', filename);
save(filename, 'L', 'n', 'r', 'GAUSSIAN_KERNEL_SIGMA', 'diffvols');
