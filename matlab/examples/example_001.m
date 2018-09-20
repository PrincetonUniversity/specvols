% EXAMPLE_001 Sample pipeline for estimating covariance from simulated data

% The number of distinct volumes. The rank of the population covariance is
% (C-1), so for higher values of C, more eigenvectors need to be estimated
% accurately to distinguish between the volumes.
C = 2;

% Image size, L-by-L. For L = 8, we can reproduce most behavior of the
% algorithm, while L = 16 is applicable to experimental data. For larger
% values of L, the memory requirements become quite significant.
L = 8;

% Number of images. The larger L is, the more images are needed to properly
% resolve the high frequencies. Runtime is linear in n.
n = 1024;

% Number of eigenvectors to estimate. This is typically (C-1), since this is
% the expected rank of a covariance matrix obtained from C distinct volumes.
% In experimental settings, however, this is of course not known in advance
% and a larger value must be selected. Here, we select num_eigs = 16 to
% get an idea of the spectrum of the estimated covariance matrix.
num_eigs = 16;

% Noise variance and spectral density. Here, we use white noise through
% `scalar_filter`, so the noise is iid with the variance specified here.
noise_var = 0;
noise_psd = scalar_filter(noise_var);

% Set of filters for the images. These objects represent the different
% filters to apply to the images during the image formation process. Here,
% we have chosen a default set of 7 standard contrast transfer functions
% (CTFs).
filters = ctf_filter();

% The ground truth volumes. Here, we generate Gaussian blobs at random
% locations as phantom volumes. These are stored as L-by-L-by-L arrays.
vols = gaussian_blob_vols(L, C);

% Parameters for the simulation.
sim_params = struct();
sim_params.n = n;
sim_params.noise_psd = noise_psd;
sim_params.filters = filters;
sim_params.vols = vols;

% Options for mean estimation. Optionally, we can set the required tolerance
% of the conjugate gradient method here, as well as the maximum number of
% iterations and other parameters.
mean_est_opt = struct();
mean_est_opt.verbose = 1;

% Options for the covariance estimation. As for the mean estimation, other
% parameters may also be set here.
covar_est_opt = struct();
covar_est_opt.verbose = 1;

% Create simulation object from parameters. This doesn't actually generate
% any images, but it sets up all the parameters of the model, including
% viewing angles (rotations), translations, amplitudes, filter assignments
% and volume assignments (that is, which volume corresponds to which image).
sim = create_sim(sim_params);

% Convert simulation object to source object. All the functions below act
% on source objects, which can come from a simulation or from experimental
% data. A source is a collection of images along with some metadata, such
% as filters. In this case, the rotations, translations, and amplitudes
% are also transferred to the source object.
src = sim_to_src(sim);

% Create a basis. Here, we use a standard spherical Fourier-Bessel basis to
% represent the volumes. This leads to a better-conditioned inverse problem
% and also increases the signal-to-noise ratio (SNR) by excluding information
% outside the unit ball in the spatial domain.
basis = fb_basis(src.L*ones(1, 3));

% Precompute the basis vectors to speed up calculations. This step is optional.
basis = precompute_basis(basis);

% Estimate the noise variance. This is needed for the covariance estimation
% step below.
noise_var_est = estimate_noise_power(src);

% Estimate the mean. This uses conjugate gradient on the normal equations for
% the least-squares estimator of the mean volume. The mean volume is
% represented internally using the basis object, but the output is in the form
% of an L-by-L-by-L array.
mean_est = estimate_mean(src, basis, mean_est_opt);

% Estimate the covariance. Again, this uses conjugate gradient on normal
% equations, but here on the least-squares estimator for the covariance.
% As before, the covariance matrix is represented in the specified basis
% during estimation, but the output is in the form of a volume matrix, that is
% an L-by-L-by-L-by-L-by-L-by-L array.
covar_est = estimate_covar(src, mean_est, noise_var_est, basis, covar_est_opt);

% Extract the top eigenvectors and eigenvalues of the covariance estimate.
% Since we know the population covariance is low-rank, we are only interested
% in the top eigenvectors.
[eigs_est, lambdas_est] = mdim_eigs(covar_est, num_eigs, 'la');

% Truncate the eigendecomposition. Since we know the true rank of the
% covariance matrix, we enforce it here.
eigs_est_trunc = eigs_est(:,:,:,1:C-1);
lambdas_est_trunc = lambdas_est(1:C-1,1:C-1);

% Estimate the coordinates in the eigenbasis. Given the images, we find the
% coordinates in the basis that minimize the mean squared error, given the
% (estimated) covariances of the volumes and the noise process.
coords_est = src_wiener_coords(src, mean_est, eigs_est_trunc, ...
    lambdas_est_trunc, noise_var_est);

% Cluster the coordinates using k-means. Again, we know how many volumes we
% expect, so we can use this parameter here. Typically, one would take the
% number of clusters to be one plus the number of eigenvectors extracted.
[centers, idx] = cluster_kmeans(coords_est, size(sim.vols, 4));

% Assign the cluster centroids to the different images. Since we expect a
% discrete distribution of volumes (and therefore of coordinates), we assign
% the centroid coordinate to each image that belongs to that cluster.
clustered_coords_est = centers(:,idx);

% Evaluate performance of mean estimation.
mean_perf = sim_eval_mean(sim, mean_est);

% Evaluate performance of covariance estimation. We also evaluate the truncated
% eigendecomposition. This is expected to be a closer approximation since it
% imposes an additional low-rank condition on the estimate.
covar_perf = sim_eval_covar(sim, covar_est);
eigs_perf = sim_eval_eigs(sim, eigs_est_trunc, lambdas_est_trunc);

% Evaluate clustering performance. We also evaluate the coordinates estimated
% after the clustering.
clustering_perf = sim_eval_clustering(sim, idx);
coords_perf = sim_eval_coords(sim, mean_est, eigs_est_trunc, clustered_coords_est);

% Output population covariance spectrum.
[~, lambdas_true] = sim_eigs(sim);
fprintf('Population covariance spectrum:\n');
for k = 1:size(lambdas_true, 1)
    fprintf('    %8g\n', lambdas_true(k,k));
end
fprintf('\n');

% Output estimated covariance spectrum.
fprintf('Population covariance spectrum:\n');
for k = 1:num_eigs
    fprintf('    %8g\n', lambdas_est(k,k));
end
fprintf('\n');

% Output performance results.
fmt_str = '%-40s%8g\n';
fprintf(fmt_str, 'Mean (rel. error):', mean_perf.rel_err);
fprintf(fmt_str, 'Mean (correlation):', mean_perf.corr);
fprintf(fmt_str, 'Covariance (rel. error):', covar_perf.rel_err);
fprintf(fmt_str, 'Covariance (correlation):', covar_perf.corr);
fprintf(fmt_str, 'Eigendecomposition (rel. error):', eigs_perf.rel_err);
fprintf(fmt_str, 'Eigendecomposition (cos theta_max):', eigs_perf.cos_principal_angles(C-1));
fprintf(fmt_str, 'Clustering (accuracy):', clustering_perf.accuracy);
fprintf(fmt_str, 'Coordinates (mean rel. error):', mean(coords_perf.rel_err));
fprintf(fmt_str, 'Coordinates (mean correlation):', mean(coords_perf.corr));
fprintf('\n');
