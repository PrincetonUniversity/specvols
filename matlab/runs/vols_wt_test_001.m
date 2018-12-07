%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tests some of the basic functionality of the `vols_wt`
% estimation. Ideally, all of these tests should yield errors close to machine
% precision when rel_tolerance is adjusted to lower values.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries: Set up the simulations, estimation options, etc.

% Set some parameters for the simulation.

n = 1024;
L = 16;

% Control the datatype of the simulation and what tolerance we should have CG
% converge to.

precision = 'single';
rel_tolerance = 1e-3;

% Initialize basis and options structs.

basis = ffb_basis(L*ones(1, 3));

sim_params = struct();
sim_params.n = n;
sim_params.L = L;
sim_params.noise_psd = scalar_filter(0);

mean_est_opt = struct();
mean_est_opt.rel_tolerance = rel_tolerance;
mean_est_opt.max_iter = 250;
mean_est_opt.precision = precision;

vols_wt_est_opt = struct();
vols_wt_est_opt.rel_tolerance = rel_tolerance;
vols_wt_est_opt.max_iter = 250;
vols_wt_est_opt.precision = precision;

% Define output function and print some of the parameters.

print_value = @(label, err)(fprintf('%-30s%.8e\n', [label ':'], err));

print_value('Number of images (n)', n);
print_value('Image size (L)', L);
print_value('Accuracy (rel_tolerance)', rel_tolerance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 1: Fit the mean volume in a single-volume simulation (C = 1).

% Set volumes, create simulation, and convert to source.

vols = gaussian_blob_vols(L, 1, [], [], precision);
vols = basis_project(basis, vols);
sim_params.vols = vols;

sim = create_sim(sim_params);
src = sim_to_src(sim);

% In the C = 1 case, we expect to recover the mean volume exactly (provided n
% is large enough), so let's get the true mean from the simulation.

mean_true = sim_mean(sim);

% All the weights are equal here to get the mean volume.

wts = ones(1, n)/sqrt(n);

% Estimate the volume.

vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

% Again, this should be close to the relative tolerance specified times the
% condition number of the forward operator (which is on the order of 10).

print_value('Test 1', anorm(mean_true-vols_wt_est)/anorm(mean_true))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 2: Fit the mean volume in a five-volume simulation (C = 5).

% Set volumes, etc. (see above).

vols = basis_project(basis, gaussian_blob_vols(L, 5, [], [], precision));
sim_params.vols = vols;

sim = create_sim(sim_params);
src = sim_to_src(sim);

% Again, weights are all equal.

wts = ones(1, n)/sqrt(n);

% Since C > 1, we don't expect the estimate to coincide with the true mean
% volume. However, the case where all weights are equal should yield the same
% estimate as the tradtional mean estimator (`estimate_mean`), so let's
% compute that.

mean_est = estimate_mean(src, basis, mean_est_opt);

% Again, estimate.

vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

% This should be close to the relative tolerance times the condition number.

print_value('Test 2', anorm(mean_est-vols_wt_est)/anorm(mean_est))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 3: For a zero-mean volume configuration (i.e., `sim_mean` is the zero
% volume) with two volumes (C = 2), recover the volumes.

% Set volumes, etc. Here, we do some preprocessing to ensure that
% the mean of the volumes is zero.

vols = gaussian_blob_vols(L, 2, [], [], precision);
vols = bsxfun(@minus, vols, mean(vols, 4));
vols = basis_project(basis, vols);
sim_params.vols = vols;

sim = create_sim(sim_params);
src = sim_to_src(sim);

% Set the weights so that volume 1 gets assigned weight 1 and volume 2 gets
% assigned weight -1.

wts = -2*(sim.states-3/2)/sqrt(n);

% Estimate.

vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

% Since we set volume 1 to get weight 1, our estimate should correspond to
% that one (volume 2 is just volume 1 multiplied by -1, so we're good).

print_value('Test 3a', anorm(vols(:,:,:,1)-vols_wt_est)/anorm(vols(:,:,:,1)))

% Let's add a "dummy" mean volume that we're estimating to see if that changes
% the behavior. Since the mean volume is zero here, we should just get zero
% for the first estimated volume, while the second gets us what we had before.
% First weight is thus 1 for all images, while second varies with the class,
% as before.

wts = [ones(1, n); -2*(sim.states-3/2)]/sqrt(n);

% Estimate.

vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

% Since we have zero mean, the norm of the first estimated volume should be
% small.

print_value('Test 3b', anorm(vols_wt_est(:,:,:,1)))

% Again, since the second estimated volume corresponds to the "difference"
% from the mean, this should be close to volume 1 in the simulation.
print_value('Test 3c', anorm(vols(:,:,:,1)-vols_wt_est(:,:,:,2))/anorm(vols(:,:,:,1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 4: A standard two-volume configuration (C = 2).

% Set the volumes, etc. Nothing special here.

vols = gaussian_blob_vols(L, 2, [], [], precision);
vols = basis_project(basis, vols);
sim_params.vols = vols;

sim = create_sim(sim_params);
src = sim_to_src(sim);

% As in the last part of the previous test, we set the first weight to one for
% all the images, while the second gives us the "coordinates" in terms of the
% difference volume.

wts = [ones(1, n); -2*(sim.states-3/2)]/sqrt(n);

% As before, we should be able to get the true mean and the true difference
% volume here, so let's compute them from the simulation.

mean_true = sim_mean(sim);
diff_true = vols(:,:,:,1)-mean_true;

vols_wt_true = cat(4, mean_true, diff_true);

% While we're at it, let's check whether this set of true volumes satisfies
% the normal equations. First, we calculate the kermat and the right-hand
% side.

kermat_f = src_vols_wt_kermat(src, wts*sqrt(n), vols_wt_est_opt);
b_coeff = src_vols_wt_backward(src, basis, wts*sqrt(n), vols_wt_est_opt);

% Calculate the basis decomposition of these true volumes and apply the
% kernel to those coefficients.

vols_wt_true_coeff = basis_expand(basis, vols_wt_true);

b_true_coeff = apply_vols_wt_kernel(vols_wt_true_coeff, kermat_f, basis, ...
    vols_wt_est_opt);

% This should be (really) small (i.e., close to eps(precision)), since the
% true volumes should satisfy the normal equations exactly.

print_value('Test 4a', anorm(b_coeff-b_true_coeff)/anorm(b_coeff));

% Now try estimating the volumes from the right-hand side we computed from the
% images.

vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

% Since the true volumes satisfy the normal equations exactly, we expect that,
% by inverting those normal equations, we should get the same volumes, up to
% the relative tolerance times the condition number.

print_value('Test 4b', anorm(vols_wt_true-vols_wt_est)/anorm(vols_wt_true));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 5: Now do five volumes (C = 5).

% Set the volumes, etc. Nothing to see.

vols = gaussian_blob_vols(L, 5, [], [], precision);
vols = basis_project(basis, vols);
sim_params.vols = vols;

sim = create_sim(sim_params);
src = sim_to_src(sim);

% Things get a little more complicated for C > 2, so let's just use the
% built-in mean and eigenvector functionality and get the coordinates in that
% (affine) space.

mean_true = sim_mean(sim);
eigs_true = sim_eigs(sim);

vol_coords = sim_vol_coords(sim, mean_true, eigs_true);

% Like before, the first weight is constant and should recover the mean
% volume. The remaining C-1 weights correspond to the volume coordinates
% obtained above.

wts = [ones(1, n); vol_coords(:,sim.states)]/sqrt(n);

% This next step is not strictly necessary. Everything is well defined without
% it. However, if we do not normalize the weights to have the same norm, the
% forward operator can easily become ill-conditioned. It will still yield the
% correct solution, but it will take a lot more iterations to converge.
% Fortunately, normalizing the weights is equivalent to multiplying the
% volumes by the same factor, so let's do that.

nrm = sqrt(sum(abs(wts).^2, 2));
wts = bsxfun(@times, wts, 1./nrm);

eigs_true_nrm = bsxfun(@times, eigs_true, permute(nrm(2:end), [2 3 4 1]));

% Combine the mean and ("normalized") eigenvectors to get the set of "true"
% volumes.

vols_wt_true = cat(4, mean_true, eigs_true_nrm);

% Let's check the normal equations again.

kermat_f = src_vols_wt_kermat(src, wts*sqrt(n), vols_wt_est_opt);
b_coeff = src_vols_wt_backward(src, basis, wts*sqrt(n), vols_wt_est_opt);

vols_wt_true_coeff = basis_expand(basis, vols_wt_true);

b_true_coeff = apply_vols_wt_kernel(vols_wt_true_coeff, kermat_f, basis, ...
    vols_wt_est_opt);

% As before, this should be really tiny.

print_value('Test 5a', anorm(b_coeff-b_true_coeff)/anorm(b_coeff));

% Finally, we estimate the volumes from the normal equations.

vols_wt_est = estimate_vols_wt(src, basis, wts, vols_wt_est_opt);

% The estimated volumes should equal the true volumes, up to the relative
% tolerance times the condition number. Note that, even with the
% normalization, the condition number for larger C is still quite bit, so we
% don't expect to hit the same accuracy as for the previous tests.

print_value('Test 5b', anorm(vols_wt_true-vols_wt_est)/anorm(vols_wt_true));
