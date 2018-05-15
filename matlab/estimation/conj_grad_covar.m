% CONJ_GRAD_COVAR Solve for volume covariance using conjugate gradient method
%
% Usage
%    [covar_est, cg_info] = conj_grad_covar(covar_kernel_f, covar_b, ...
%        basis, precond_kernel_f, covar_est_opt);
%
% Input
%    covar_kernel_f: The non-centered Fourier transform of the covariance
%       projection-backprojection operator obtained from `src_covar_kernel`.
%    covar_b: An volume matrix of size L-by-L-by-L-by-L-by-L-by-L containing
%       the backprojected covariance estimates obtained from
%       `src_covar_backward`.
%    basis: A basis object used for representing the volumes (default
%       `dirac_basis(L*ones(1, 3))`).
%    precond_kernel_f: If not empty, the Fourier transform of a kernel that is
%       used to precondition the covariance projection-backprojection operator
%       (default empty).
%    covar_est_opt: An options structure. No options are used in the function
%       itself, but it is passed on to the `conj_grad` function, so any options
%       for that function should be specified here.
%
% Output
%    covar_est: A volume matrix of size L-by-L-by-L-by-L-by-L-by-L containing
%       the least-squares estimate obtained by solving the eqation L(X) = B,
%       where A is the covariance projection-backprojection operator represented
%       by `covar_kernel_f` and B are the backprojected covariance estimates in
%       `covar_b`. The equation is solved using the conjugate gradient method
%       implemented by `conj_grad`.
%    cg_info: A structure containing information about the conjugate gradient
%       method, such as residuals, objectives, etc. See the documentation of
%       `conj_grad` for more details.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [covar_est, cg_info] = conj_grad_covar(covar_kernel_f, covar_b, ...
    basis, precond_kernel_f, covar_est_opt)

    if nargin < 3
        basis = [];
    end

    if nargin < 4
        precond_kernel_f = [];
    end

    if nargin < 5 || isempty(covar_est_opt)
        covar_est_opt = struct();
    end

    L = size(covar_b, 1);

    if isempty(basis)
        basis = dirac_basis(L*ones(1, 3));
    end

    if ndims(covar_b) ~= 6 || any(size(covar_b) ~= L)
        error(['Input `covar_b` must be a volume matrix of size ' ...
            'L-by-L-by-L-by-L-by-L-by-L.']);
    end

    if ~is_basis(basis) || any(basis.sz ~= L*ones(1, 3))
        error(['Input `basis` must be a basis object representing ' ...
            'volumes of size L-by-L-by-L.']);
    end

    kernel_fun = @(covar_coeff, kernel_f)(symmat_to_vec_iso( ...
        apply_covar_kernel(vec_to_symmat_iso(covar_coeff), ...
        kernel_f, basis, covar_est_opt)));

    covar_kernel_fun = @(covar_coeff)(kernel_fun(covar_coeff, covar_kernel_f));

    if ~isempty(precond_kernel_f)
        precond_kernel_fun = ...
            @(covar_coeff)(kernel_fun(covar_coeff, precond_kernel_f));

        covar_est_opt.preconditioner = precond_kernel_fun;
    end

    covar_b_coeff = symmat_to_vec_iso(basis_mat_evaluate_t(basis, covar_b));

    [covar_est_coeff, ~, cg_info] = ...
        conj_grad(covar_kernel_fun, covar_b_coeff, covar_est_opt);

    covar_est = basis_mat_evaluate(basis, vec_to_symmat_iso(covar_est_coeff));
end
