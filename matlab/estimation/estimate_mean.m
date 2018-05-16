% ESTIMATE_MEAN Estimate mean using least squares and conjugate gradients
%
% Usage
%    [mean_est, cg_info] = estimate_mean(src, basis, mean_est_opt);
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    basis: A basis object used for representing the volumes (default
%       dirac_basis(L*ones(1, 3))).
%    mean_est_opt: A struct containing the fields:
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
%    mean_est: The estimated mean volume, in the form of an L-by-L-by-L array.
%       It minimizes the objective
%
%          1/n sum_{s=1}^n |P_s x - y_s|^2 + lambda | x |^2,
%
%       where x is the volume, P_s are the imaging mappings (basis evaluation,
%       projection, CTF filtering, translation, scaling), lambdas is the
%       regularization parameter `mean_est_opt.regularizer`, and y_s are the
%       observed images. This is achieved by forming the normal equations and
%       solving them using the preconditioned conjugate gradient method.
%    cg_info: A structure containing information about the conjugate gradient
%       method, such as residuals, objectives, etc. See the documentation of
%       `conj_grad` for more details.
%
% See also
%    src_mean_kernel, src_mean_backward

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [mean_est, cg_info] = estimate_mean(src, basis, mean_est_opt)
    if nargin < 2
        basis = [];
    end

    if nargin < 3 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    L = src.L;
    n = src.n;

    mean_est_opt = fill_struct(mean_est_opt, ...
        'preconditioner', 'circulant', ...
        'precision', 'single');

    if isempty(basis)
        basis = dirac_basis(L*ones(1, 3));
    end

    kernel_f = src_mean_kernel(src, mean_est_opt);

    precond_kernel_f = [];

    if ischar(mean_est_opt.preconditioner)
        if strcmp(mean_est_opt.preconditioner, 'none')
            precond_kernel_f = [];
        elseif strcmp(mean_est_opt.preconditioner, 'circulant')
            precond_kernel_f = circularize_kernel(kernel_f);
        else
            error('Invalid preconditioner type.');
        end

        % Reset so this is not used by the `conj_grad` function.
        mean_est_opt.preconditioner = @(x)(x);
    end

    mean_b = src_mean_backward(src, mean_est_opt);

    [mean_est, cg_info] = conj_grad_mean(kernel_f, mean_b, basis, ...
        precond_kernel_f, mean_est_opt);
end
