% Script to estimate the condition number `kappa` of the mean projection-
% backprojection operator.
%
% It calculates `kappa` first using the iterative method, which feeds the
% convolution into `eigs` to get the largest eigenvalue, then defines its
% inverse using `conj_grad` which is then fed into `eigs` to get the smallest
% eigenvalue. The true value of `kappa` is also calculated by forming the full
% matrix corresponding to the convolution and using `cond`.
%
% Note that the accuracy of `lambda_min` (and consequently, that of
% `kappa_eigs`) depends on the convergence of the CG. Specifically, if we set
% `mean_est_opt.rel_tolerance` to a lower value, we will get a more accurate
% estimate for `lambda_min`. On the other hand, this will take a longer time
% to converge, especially if the operator is ill-conditioned. In other words,
% don't try to estimate condition numbers of ill-conditioned operators. That
% estimation problem is itself ill-conditioned.

% sim = create_sim();
% src = sim_to_src(sim);

% mean_kernel_f = src_mean_kernel(src);


% basis = fb_basis(src.L*ones(1, 3));

mean_est_opt = struct();
mean_est_opt.max_iter = 10;
mean_est_opt.rel_tolerance = 1e-3;

eigs_opt = struct();
eigs_opt.isreal = true;
eigs_opt.issym = true;

% A_fun = @(x)(apply_mean_kernel(x, mean_kernel_f, basis, mean_est_opt));
A_fun = @(x)(apply_vols_wt_kernel(x,kermat_f,basis,mean_est_opt));

Ainv_fun = @(x)(conj_grad_vols_wt(A_fun, x, mean_est_opt));

lambda_max = eigs(A_fun, basis.count, 1, 'lm', eigs_opt);
lambda_min = eigs(Ainv_fun, basis.count, 1, 'sm', eigs_opt);

kappa_eigs = lambda_max/lambda_min;

% A = kernel_to_toeplitz(mean_kernel_f);
l = size(kermat_f,1);
r = size(kermat,4);
A = zeros([(size(kermat_f,1)^3*size(kermat,4)) (size(kermat_f,1)^3*size(kermat,5))]);
for i = 1:size(kermat_f,4)
    for j = 1:size(kermat_f,5
        A() = basis_mat_evaluate_t(basis, A);
    end
end

kappa_true = cond(A);

fprintf('Estimated condition number:   %g\n', kappa_eigs);
fprintf('True condition number:        %g\n', kappa_true);