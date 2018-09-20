% ESTIMATE_COVAR Estimate covariance using least squares
%
% Usage
%    [covar_est, cg_info] = estimate_covar(src, mean_vol, noise_var, ...
%        basis, covar_est_opt);
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    mean_vol: The (estimated) mean volume of the source. This can be estimated
%       using `estimate_mean`.
%    noise_var: The variance of the noise.
%    basis: A basis object used for representing the volumes (default
%       dirac_basis(L*ones(1, 3))).
%    covar_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double' or
%             'single' (default).
%          - 'batch_size': The size of the batches in which to compute the
%             kernel (default 512).
%          - 'regularizer': The regularizer parameter for the least-squares
%             problem. This is a positive number that determines how much
%             the least-squares solution is to be regularized (default 0).
%          - 'preconditioner': One of the following values specifying the
%             preconditioner for the conjugate gradient method:
%                - 'none': No preconditioner is used.
%                - 'circulant': Uses `circularize_kernel` to obtain a
%                   circulant approximation to the projection-backprojection
%                   kernel whose inverse is then used as a preconditioner
%                   (default).
%                - function handle: In this case, the function handle is passed
%                   on directly to the `conj_grad` function.
%       The struct is also passed on to the `conj_grad` function, so any
%       options to that function should be passed here.
%
% Output
%    covar_est: The estimated volume covariance matrix, in the form of an
%       L-by-L-by-L-by-L-by-L-by-L array. It minimizes the objective
%
%          1/n sum_{s=1}^n | P_s X P_s^T + noise_var I - 
%             (y_s - P_s mean_vol) (y_s - P_s mean_vol)^T |_F^2 +
%             lambda * | X |_F^2,
%
%       where X is the covariance matrix estimate, P_s are the imaging
%       mappings (basis evaluation, projection, CTF filtering, translation,
%       scaling), lambda is the regularization parameter
%       `covar_est_opt.regularizer`, and y_s are the observed images. The
%       objective is minimized by forming the normal equations and solving
%       them using the preconditioned conjugate gradient method.
%    cg_info: A structure containing information about the conjugate gradient
%       method, such as residuals, objectives, etc. See the documentation of
%       `conj_grad` for more details.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function covar_est = estimate_covar(src, mean_vol, noise_var, basis, ...
    covar_est_opt)

    if nargin < 4
        basis = [];
    end

    if nargin < 5 || isempty(covar_est_opt)
        covar_est_opt = struct();
    end

    L = src.L;
    n = src.n;

    covar_est_opt = fill_struct(covar_est_opt, ...
        'preconditioner', 'circulant', ...
        'precision', 'single');

    if isempty(basis)
        basis = dirac_basis(L*ones(1, 3));
    end

    kernel_f = src_covar_kernel(src, covar_est_opt);

    precond_kernel_f = [];

    if ischar(covar_est_opt.preconditioner)
        if strcmp(covar_est_opt.preconditioner, 'none')
            precond_kernel_f = [];
        elseif strcmp(covar_est_opt.preconditioner, 'circulant')
            precond_kernel_f = circularize_kernel(kernel_f);
        else
            error('Invalid preconditioner type.');
        end

        % Reset so this is not used by the `conj_grad` function.
        covar_est_opt.preconditioner = @(x)(x);
    end

    b_coeff = src_covar_backward(src, basis, mean_vol, noise_var, ...
       covar_est_opt);

    [covar_est_coeff, cg_info] = conj_grad_covar(kernel_f, b_coeff, basis, ...
        precond_kernel_f, covar_est_opt);

    covar_est = basis_mat_evaluate(basis, covar_est_coeff);

    covar_est = vecmat_to_volmat(make_symmat(volmat_to_vecmat(covar_est)));
end
