%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From altmin_test_001.m
% Instead of perturbing the true weights, we initialize with completely random
% weights (`init_type = 1`). For comparison, we have the results when
% initializing with true weights (`init_type = 0`). The expectation is that
% this will converge to something like that true (global) least-squares
% solution. With that in mind, the randomly initialized results are also
% compared with the final result of the true initialization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('n'), n = 1024; end
if ~exist('L'), L = 16; end

% Control the datatype of the simulation and what tolerance we should have CG
% converge to.

precision = 'single';
rel_tolerance = 1e-6;

% Set noise on the images.

if ~exist('im_noise_var'), im_noise_var = 3.36e-2; end

% Set the number of alternating minimization iterations.

altmin_iter = 50;

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
vols_wt_est_opt.batch_size = 2^13;

% Define some functions used later for output.

cos_max_principal_angle = @(vols1, vols2)( ...
    min(cos_principal_angles(vol_to_vec(vols1), vol_to_vec(vols2))));
print_value = @(label, err)(fprintf('%-40s%.8e\n', [label ':'], err));

% Print basic parameters of the experiment.

print_value('Number of images (n)', n);
print_value('Image size (L)', L);
print_value('Accuracy (rel_tolerance)', rel_tolerance);
print_value('Image noise variance', im_noise_var);

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

% Run alternating minimization starting from "true" weights and random weights.

for init_type = [0 1]
    if init_type == 0
        fprintf('Initializing with true weights.\n');
        wts = wts_true;
    elseif init_type == 1
        fprintf('Initializing with random weights.\n');
        wts = wts_true;
        wts(2:end,:) = randn(size(wts_true,1)-1, size(wts, 2));
    end

    nrm = sqrt(sum(abs(wts).^2, 2));
    wts = bsxfun(@times, wts, 1./nrm);

    vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

    corr(1,init_type+1) = ...
        cos_max_principal_angle(vols_wt_true, vols_wt_est);

    print_value('Correlation (init.)', corr(1,init_type+1));

    if init_type == 1
        corr(1,3) = ...
            cos_max_principal_angle(vols_wt_est_itrue, vols_wt_est);
    end

    for k = 1:altmin_iter
        [Q, ~] = qr(vol_to_vec(vols_wt_est(:,:,:,2:end)), 0);
        vols_wt_est(:,:,:,2:end) = vec_to_vol(Q);

        coords = src_wiener_coords(src, ...
            vols_wt_est(:,:,:,1), vols_wt_est(:,:,:,2:end));

        wts(2:end,:) = coords;

        nrm = sqrt(sum(abs(wts).^2, 2));
        wts = bsxfun(@times, wts, 1./nrm);

        vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

        corr(k+1,init_type+1) = ...
            cos_max_principal_angle(vols_wt_true, vols_wt_est);

        print_value('Correlation (iteration)', corr(k+1,init_type+1));

        if init_type == 1
            corr(k+1,3) = ...
                cos_max_principal_angle(vols_wt_est_itrue, vols_wt_est);
        end
    end

    if init_type == 0
        vols_wt_est_itrue = vols_wt_est;
    elseif init_type == 1
        vols_wt_est_irand = vols_wt_est;
    end
end

print_value('Correlation true vs. rand init.', corr(end,3));
