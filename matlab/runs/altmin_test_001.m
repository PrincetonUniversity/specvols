%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is a proof-of-concept for the alternating minimization algorithm.
% It sets up a heterogeneity problem, finds the true weights, and attempts to
% find the correct volume space by alternating minimization on that space and
% the weights.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1024;
L = 16;

% Control the datatype of the simulation and what tolerance we should have CG
% converge to.

precision = 'single';
rel_tolerance = 1e-6;

% Set noise on the images and "noise" in the initial coordinates.

im_noise_var = 3.36e-2;
coord_noise_var = 100;

% Set the number of alternating minimization iterations.

altmin_iter = 20;

% Initialize basis and options structs.

basis = ffb_basis(L*ones(1, 3));

sim_params = struct();
sim_params.n = n;
sim_params.L = L;
sim_params.noise_psd = scalar_filter(im_noise_var);

vols = gaussian_blob_vols(L, 3, [], [], precision);
vols = basis_project(basis, vols);
sim_params.vols = vols;

vols_wt_est_opt = struct();
vols_wt_est_opt.rel_tolerance = rel_tolerance;
vols_wt_est_opt.max_iter = 250;
vols_wt_est_opt.precision = precision;

% Define some functions used later for output.

cos_max_principal_angle = @(vols1, vols2)( ...
    min(cos_principal_angles(vol_to_vec(vols1), vol_to_vec(vols2))));
print_value = @(label, err)(fprintf('%-30s%.8e\n', [label ':'], err));

% Print basic parameters of the experiment.

print_value('Number of images (n)', n);
print_value('Image size (L)', L);
print_value('Accuracy (rel_tolerance)', rel_tolerance);
print_value('Image noise variance', im_noise_var);
print_value('Coordinate noise variance', coord_noise_var);

% Create simulation, source, and cache it. Since we're going to do a lot of
% backprojections with different weights, we save time by doing this here.

sim = create_sim(sim_params);
src = sim_to_src(sim);

src = cache_src(src);

% Get true mean and eigenvectors to generate "true" weights and volume space.

mean_true = sim_mean(sim);
eigs_true = sim_eigs(sim);

vol_coords = sim_vol_coords(sim, mean_true, eigs_true);

wts_true = [ones(1, n); vol_coords(:,sim.states)]/sqrt(n);

vols_wt_true = cat(4, mean_true, eigs_true);

% Normalize weights to improve convergence.

nrm = sqrt(sum(abs(wts_true).^2, 2));
wts_true = bsxfun(@times, wts_true, 1./nrm);

% Calculate volumes that correspond to the "true" weights.

vols_wt_est = estimate_vols_wt(src, basis, wts_true, vols_wt_est_opt);

print_value('Correlation (true weights)', ...
    cos_max_principal_angle(vols_wt_true, vols_wt_est));

% Perturb the "true" weights to get initial weights for the algorithm.

wts = wts_true;
wts(2:end,:) = wts(2:end,:) + ...
    sqrt(coord_noise_var)*randn(size(wts, 1)-1, size(wts, 2));

% Normalize weights and calculate volumes corresponding to the initial weights.

nrm = sqrt(sum(abs(wts).^2, 2));
wts = bsxfun(@times, wts, 1./nrm);

vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

print_value('Correlation (init.)', ...
    cos_max_principal_angle(vols_wt_true, vols_wt_est));

% Alternate between calculating weights and estimating volumes.

for k = 1:altmin_iter
    coords = src_wiener_coords(src, ...
        vols_wt_est(:,:,:,1), vols_wt_est(:,:,:,2:end));

    wts(2:end,:) = coords;

    nrm = sqrt(sum(abs(wts).^2, 2));
    wts = bsxfun(@times, wts, 1./nrm);

    vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

    print_value('Correlation (iteration)', ...
        cos_max_principal_angle(vols_wt_true, vols_wt_est));
end
