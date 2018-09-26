% ESTIMATE_VOLS_WT Estimate vols for weighted means using least squares and
% conjugate gradients
%
% Usage
%    [vols_wt_est, cg_info] = estimate_vols_wt(src, basis, wts, vols_wt_est_opt)
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    basis: A basis object used for representing the volumes (default
%       dirac_basis(L*ones(1, 3))).
%    wts: a matrix whose columns contain the imaging weightings for each of
%       the "diffusion volumes" (we need a better, more general way of
%       saying this).
%
%    vols_wt_est_opt: A struct containing the fields:
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
%    vols_wt_est: The estimated mean volume, in the form of an L-by-L-by-L array.
%       It minimizes the objective
%
%          1/n sum_{s=1}^n |\P_s (sum_{l=0}^{r-1} a_l)|^2 +
%               sum_{l=0}^{r-1} lambda_l |a_l|^2
%
%       where x is the volume, P_s are the imaging mappings (basis evaluation,
%       projection, CTF filtering, translation, scaling), lambdas is the
%       regularization parameter `vols_wt_est_opt.regularizer`, and y_s are the
%       observed images. This is achieved by forming the normal equations and
%       solving them using the preconditioned conjugate gradient method.
%    cg_info: A structure containing information about the conjugate gradient
%       method, such as residuals, objectives, etc. See the documentation of
%       `conj_grad` for more details.
%
% See also
%    src_vols_wt_kernel, src_vols_wt_backward

% Author
%    Joakim Anden <janden@flatironinstitute.org>
%    Amit Halevi <ahalevi@princeton.edu>

function [vols_wt_est, cg_info] = estimate_vols_wt(src, basis, wts, vols_wt_est_opt)
    if nargin < 2
        basis = [];
    end

    if nargin < 4 || isempty(vols_wt_est_opt)
        vols_wt_est_opt = struct();
    end

    L = src.L;
    n = src.n;

%     vols_wt_est_opt = fill_struct(vols_wt_est_opt, ...
%         'preconditioner', 'circulant', ...
%         'precision', 'single');
    vols_wt_est_opt = fill_struct(vols_wt_est_opt, ...
        'preconditioner', 'none', ...
        'precision', 'single');
    if isempty(basis)
        basis = dirac_basis(L*ones(1, 3));
    end

    kermat_f = sqrt(n^2) * src_vols_wt_kermat_slices(src, wts, vols_wt_est_opt);

    precond_kermat_f = [];

    if ischar(vols_wt_est_opt.preconditioner)
        if strcmp(vols_wt_est_opt.preconditioner, 'none')
            precond_kermat_f = [];
        elseif strcmp(vols_wt_est_opt.preconditioner, 'circulant')
            precond_kermat_f = circularize_kernel(kermat_f, 3);
        else
            error('Invalid preconditioner type.');
        end

        % Reset so this is not used by the `conj_grad` function.
        vols_wt_est_opt.preconditioner = @(x)(x);
    end

    vols_wt_b_coeff = sqrt(n) * src_vols_wt_backward_slices(src, basis, wts, vols_wt_est_opt);

    
    [vols_wt_est_coeff, cg_info] = conj_grad_vols_wt(kermat_f, vols_wt_b_coeff, ...
        basis, precond_kermat_f, vols_wt_est_opt);

    vols_wt_est = basis_evaluate(basis, vols_wt_est_coeff);
end
