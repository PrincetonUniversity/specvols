% Script to estimate the condition number `kappa` of the preconditioned 
% mean projection-backprojection operator.
%
% It calculates `kappa` first using the iterative method, which defines the
% operator consisting of convolution by the predonditioner followed by
% convolution with the projection-backprojection kernel. The inverse of this
% operator is also defined using `conj_grad`. These operators are then fed
% into `fun_svds`, which calculates the largest singular values of an operator
% given as a function handle. Since the smallest singular value of our given
% operator is the reciprocal of the largest singular value of its inverse,
% we get both of these by applying `fun_svds` to the operator and its inverse.
% The true value of `kappa` is also calculated by forming the full
% matrix corresponding to the convolutions, multiplying them, and using
% `cond`.
%
% As before, the accuracy of `sigma_min` (and consequently, that of
% `kappa_svds`) depends on the convergence of the CG. Specifically, if we set
% `mean_est_opt.rel_tolerance` to a lower value, we will get a more accurate
% estimate for `sigma_min`.

% sim = create_sim();
% src = sim_to_src(sim);
% 
% mean_kernel_f = src_mean_kernel(src);
% 
% precond_kernel_f = circularize_kernel(mean_kernel_f);
% 
% basis = fb_basis(src.L*ones(1, 3));
% 
% mean_est_opt = struct();
% mean_est_opt.max_iter = 1000;
% mean_est_opt.rel_tolerance = 1e-9;

eigs_opt = struct();
eigs_opt.isreal = true;
eigs_opt.issym = true;

wts = dmap_coords;
kermat_f = src_vols_wt_kermat(src, wts, vols_wt_est_opt);

A_fun = @(x)(apply_vols_wt_kernel(x, kermat_f, basis, vols_wt_est_opt));
% C_fun = @(x)(apply_vols_wt_kernel(x, 1./precond_kermat_f, basis, vols_wt_est_opt));
C_fun = @(x)(x)

Ainv_fun = @(x)(conj_grad(A_fun,x, vols_wt_est_opt));
Cinv_fun = @(x)(conj_grad(C_fun,x, vols_wt_est_opt));

operator_fun = @(x)(C_fun(A_fun(x)));
operator_t_fun = @(x)(A_fun(C_fun(x)));

inverse_fun = @(x)(Ainv_fun(Cinv_fun(x)));
inverse_t_fun = @(x)(Cinv_fun(Ainv_fun(x)));

[~, operator_sigma_max] = fun_svds(operator_fun, operator_t_fun, ...
    basis.count*num_dmap_coords, basis.count*num_dmap_coords, 1, eigs_opt);
[~, inverse_sigma_max] = fun_svds(inverse_fun, inverse_t_fun, ...
    basis.count*num_dmap_coords, basis.count*num_dmap_coords, 1, eigs_opt);

sigma_max = operator_sigma_max;
sigma_min = 1./inverse_sigma_max;

toc

kappa_svds = sigma_max/sigma_min

% A = kernel_to_toeplitz(mean_kernel_f);
% A = basis_mat_evaluate_t(basis, A);
% 
% C = kernel_to_circulant(1./precond_kernel_f);
% C = basis_mat_evaluate_t(basis, C);
% 
% kappa_true = cond(C*A);
% 
% fprintf('Estimated condition number:   %g\n', kappa_svds);
% fprintf('True condition number:        %g\n', kappa_true);