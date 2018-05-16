% CONJ_GRAD_MEAN Solve for mean volume using conjugate gradient method
%
% Usage
%    mean_est = conj_grad_mean(mean_kernel_f, mean_b, basis, ...
%       precond_kernel_f, mean_est_opt);
%
% Input
%    mean_kernel_f: The non-centered Fourier transform of the projection-
%       backprojection operator obtained from `src_mean_kernel`.
%    mean_b: An array of size L-by-L-by-L containing the backprojected images
%       obtained from `src_mean_backward`.
%    basis: A basis object used for representing the volumes (default
%       `dirac_basis(L*ones(1, 3))`).
%    precond_kernel_f: If not empty, the Fourier transform of a kernel that is
%       used to precondition the projection-backprojection operator (default
%       empty).
%    mean_est_opt: An options structure. No options are used in the function
%       itself, but it is passed on to the `conj_grad` function, so any options
%       for that function should be specified here.
%
% Output
%    mean_est: An array of size L-by-L-by-L containing the least-squares
%       estimate obtained by solving the equation A*x = b, where A is the
%       linear mapping represented by convolving with `mean_kernel_f` and b
%       are the backprojected images `mean_b`. The equation is solved using the
%       conjugate gradient method implemented by `conj_grad`.
%    cg_info: A structure containing information about the conjugate gradient
%       method, such as residuals, objectives, etc. See the documentation of
%       `conj_grad` for more details.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [mean_est, cg_info] = conj_grad_mean(mean_kernel_f, mean_b, basis, ...
    precond_kernel_f, mean_est_opt)

    if nargin < 3
        basis = [];
    end

    if nargin < 4
        precond_kernel_f = [];
    end

    if nargin < 5 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    L = size(mean_b, 1);

    if isempty(basis)
        basis = dirac_basis(L*ones(1, 3));
    end

    if ndims(mean_b) ~= 3 || size(mean_b, 2) ~= L || size(mean_b, 3) ~= L
        error('Input `mean_b` must be a volume of size L-by-L-by-L.');
    end

    if ~is_basis(basis) || any(basis.sz ~= L*ones(1, 3))
        error(['Input `basis` must be a basis object representing ' ...
            'volumes of size L-by-L-by-L.']);
    end

    kernel_fun = @(vol_coeff)( ...
        apply_mean_kernel(vol_coeff, mean_kernel_f, basis, mean_est_opt));

    if ~isempty(precond_kernel_f)
        precond_fun = @(vol_coeff)( ...
            apply_mean_kernel(vol_coeff, 1./precond_kernel_f, basis, ...
            mean_est_opt));

        mean_est_opt.preconditioner = precond_fun;
    end

    mean_b_coeff = basis_evaluate_t(basis, mean_b);

    [mean_est_coeff, ~, cg_info] = ...
        conj_grad(kernel_fun, mean_b_coeff, mean_est_opt);

    mean_est = basis_evaluate(basis, mean_est_coeff);
end
