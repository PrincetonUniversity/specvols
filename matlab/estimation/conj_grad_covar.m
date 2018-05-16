% CONJ_GRAD_COVAR Solve for volume covariance using conjugate gradient method
%
% Usage
%    [covar_est_coeff, cg_info] = conj_grad_covar(covar_kernel_f, ...
%        covar_b_coeff, basis, precond_kernel_f, covar_est_opt);
%
% Input
%    covar_kernel_f: The non-centered Fourier transform of the covariance
%       projection-backprojection operator obtained from `src_covar_kernel`.
%    covar_b_coeff: A matrix of backprojected covariance estimates, expressed
%       as coefficients in a basis, obtained from `src_covar_backward`.
%    basis: A basis object used for representing the volumes.
%    precond_kernel_f: If not empty, the Fourier transform of a kernel that is
%       used to precondition the covariance projection-backprojection operator
%       (default empty).
%    covar_est_opt: A struct containing the fields:
%          - 'regularizer': The regularizer parameter for the least-squares
%             problem. This is a positive number that determines how much
%             the least-squares solution is to be regularized (default 0).
%       The struct is also passed on to the `conj_grad` function, so any
%       options to that function should be passed here.
%
% Output
%    covar_est_coeff: A matrix of size `basis.count`-by-`basis.count` contain-
%       ing the basis coefficients of the least-squares estimate obtained by
%       solving the equation L(X) + lambda X = B, where L is the covariance
%       projection-backprojection operator represented by `covar_kernel_f`,
%       lambda is the regularization parameter `covar_est_opt.lambda`, and B
%       is the sum of backprojected covariance estimates in `covar_b`, ex-
%       pressed in `basis`. The equation is solved using the preconditioned
%       conjugate gradient method implemented by `conj_grad`.
%    cg_info: A structure containing information about the conjugate gradient
%       method, such as residuals, objectives, etc. See the documentation of
%       `conj_grad` for more details.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [covar_est_coeff, cg_info] = conj_grad_covar(covar_kernel_f, ...
    covar_b_coeff, basis, precond_kernel_f, covar_est_opt)

    if nargin < 4
        precond_kernel_f = [];
    end

    if nargin < 5 || isempty(covar_est_opt)
        covar_est_opt = struct();
    end

    covar_est_opt = fill_struct(covar_est_opt, ...
        'regularizer', 0);

    if ndims(covar_b_coeff) ~= 2 || ...
        size(covar_b_coeff, 1) ~= size(covar_b_coeff, 2)

        error(['Input `covar_b_coeff` must be a matrix.']);
    end

    if ~is_basis(basis) || basis.count ~= size(covar_b_coeff, 1)
        error(['Input `basis` must be a basis object compatible with ' ...
            '`covar_b_coeff`.']);
    end

    if covar_est_opt.regularizer > 0
        covar_kernel_f = covar_kernel_f + ...
            covar_est_opt.regularizer*ones(size(covar_kernel_f));

        if ~isempty(precond_kernel_f)
            precond_kernel_f = precond_kernel_f + ...
                covar_est_opt.regularizer*ones(size(precond_kernel_f));
        end
    end

    kernel_fun = @(covar_coeff, kernel_f)(symmat_to_vec_iso( ...
        apply_covar_kernel(vec_to_symmat_iso(covar_coeff), ...
        kernel_f, basis, covar_est_opt)));

    covar_kernel_fun = @(covar_coeff)(kernel_fun(covar_coeff, covar_kernel_f));

    if ~isempty(precond_kernel_f)
        precond_kernel_fun = ...
            @(covar_coeff)(kernel_fun(covar_coeff, 1./precond_kernel_f));

        covar_est_opt.preconditioner = precond_kernel_fun;
    end

    covar_b_coeff = symmat_to_vec_iso(covar_b_coeff);

    [covar_est_coeff, ~, cg_info] = ...
        conj_grad(covar_kernel_fun, covar_b_coeff, covar_est_opt);

    covar_est_coeff = vec_to_symmat_iso(covar_est_coeff);
end
