root_dir = fullfile(pkg_root(), 'data', 'fred_ionchannel');

star_file = fullfile(root_dir, '_ct25_it050_data.star');

star_opt = struct();
star_opt.pixel_size = 1.5;

L = 16;

noise_est_opt = struct();
noise_est_opt.noise_type = 'anisotropic';

mean_est_opt = struct();
mean_est_opt.verbose = 1;
mean_est_opt.rel_tolerance = 1e-3;
mean_est_opt.batch_size = 2^13;
mean_est_opt.regularizer = 2^-16;

covar_est_opt = struct();
covar_est_opt.verbose = 1;
covar_est_opt.rel_tolerance = 1e-3;
covar_est_opt.batch_size = 2^13;
covar_est_opt.regularizer = 2^-4;

tmr = tic; fprintf('Loading STAR file...');
star = load_star(star_file);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Parsing STAR file...');
src_orig = star_to_src(star, root_dir, star_opt);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Cropping and downsampling images...');
src_cropped = src_orig;
src_ds = downsample_src(src_cropped, L);
src_ds = cache_src(src_ds);
fprintf('OK (%.2f s)\n', toc(tmr));

% TODO: This is not the right way to do it. Fix.
tmr = tic; fprintf('Symmetrizing dataset...');
rots = src_ds.params.rots;
for ell = 1:3
rot_z = angles_to_rots([ell*pi/4 0 0]');
sym_rots = zeros(3, 3, src_ds.n);
for k = 1:src_ds.n
sym_rots(:,:,k) = rot_z*src_ds.params.rots(:,:,k);
end
rots = cat(3, rots, sym_rots);
end

src_aug = struct();
src_aug.type = src_type_array();
src_aug.L = src_ds.L;
src_aug.n = 4*src_ds.n;
src_aug.precision = src_ds.precision;
src_aug.images = repmat(src_ds.images, [1 1 4]);
src_aug.params = struct();
src_aug.params.rots = rots;
src_aug.params.filters = src_ds.params.filters;
src_aug.params.filter_idx = repmat(src_ds.params.filter_idx, [1 4]);
src_aug.params.offsets = repmat(src_ds.params.offsets, [1 4]);
src_aug.params.amplitudes = repmat(src_ds.params.amplitudes, [1 4]);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Whitening images...');
noise_psd_est = estimate_noise(src_aug, noise_est_opt);
src_whiten = whiten_src(src_aug, noise_psd_est);
src_whiten = cache_src(src_whiten);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Preparing basis...');
basis = fb_basis(L*ones(1, 3), L);
basis = precompute_basis(basis);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Calculating mean kernel...');
mean_kernel_f = src_mean_kernel(src_whiten, mean_est_opt);
mean_precond_kernel_f = circularize_kernel(mean_kernel_f);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Calculating mean backprojection...');
mean_b_coeff = src_mean_backward(src_whiten, basis, mean_est_opt);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Estimating mean...\n');
[mean_est_coeff, mean_cg_info] = conj_grad_mean(mean_kernel_f, mean_b_coeff, ...
    basis, mean_precond_kernel_f, mean_est_opt);
mean_est = basis_evaluate(basis, mean_est_coeff);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Estimating noise variance...');
noise_var_est = estimate_noise_power(src_whiten);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Calculating covariance kernel...');
covar_kernel_f = src_covar_kernel(src_whiten, covar_est_opt);
covar_precond_kernel_f = circularize_kernel(covar_kernel_f);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Calculating covariance backprojection...');
covar_b_coeff = src_covar_backward(src_whiten, basis, mean_est, noise_var_est, ...
    covar_est_opt);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Estimating covariance...\n');
[covar_est_coeff, covar_cg_info] = conj_grad_covar(covar_kernel_f, ...
    covar_b_coeff, basis, covar_precond_kernel_f, covar_est_opt);
covar_est = basis_mat_evaluate(basis, covar_est_coeff);
covar_est = vecmat_to_volmat(make_symmat(volmat_to_vecmat(covar_est)));
fprintf('OK (%.2f s)\n', toc(tmr));
