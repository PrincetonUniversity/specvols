root_dir = fullfile(pkg_root(), 'data', 'yan_wu');

star_file = fullfile(root_dir, 'All_class001_r1_ct1_data.star');

star_opt = struct();
star_opt.pixel_size = 1.338;

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
covar_est_opt.regularizer = 2^-14;

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

tmr = tic; fprintf('Whitening images...');
noise_psd_est = estimate_noise(src_ds, noise_est_opt);
src_whiten = whiten_src(src_ds, noise_psd_est);
src_whiten = cache_src(src_whiten);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Preparing basis...');
basis = fb_basis(L*ones(1, 3), L);
basis = precompute_basis(basis);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Calculating mean kernel...');
mean_kernel_f = src_mean_kernel(src_whiten, mean_est_opt);
precond_kernel_f = circularize_kernel(mean_kernel_f);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Calculating mean backprojection...');
mean_b_coeff = src_mean_backward(src_whiten, basis, mean_est_opt);
fprintf('OK (%.2f s)\n', toc(tmr));

tmr = tic; fprintf('Estimating mean...\n');
[mean_est_coeff, cg_opt] = conj_grad_mean(mean_kernel_f, mean_b_coeff, ...
    basis, precond_kernel_f, mean_est_opt);
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
