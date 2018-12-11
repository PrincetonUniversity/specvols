% CONJ_GRAD_VOLS_WT Solve for a set of volumes corresponding to a set of
% weightings
%
% Usage
%    function [mean_est_coeff, cg_info] = conj_grad_vols_wt(kermat_f, ...
%    vols_wt_b_coeff, basis, precond_kermat_f, vols_wt_est_opt)
%
% Input
%    kermat_f: a "matrix" of non-centered Fourier transforms of the
%       weighted projection-backprojection operator obtained from
%       src_vols_wt_backward
%    vols_wt_b_coeff: a matrix whose columns contain the weighted means of
%       backprojected images, obtained from src_vols_wt_backward.
%    basis: The basis object used for representing the volumes.
%    precond_kermat_f: If not empty, the Fourier transform of a kernel that is
%       used to precondition the projection-backprojection operator (default
%       empty).
%    vols_wt_est_opt: The struct is also passed on to the `conj_grad`
%       function, so any options to that function should be passed here.
%
% Output
%    mean_est_vols: a matrix whose columns are of length basis.count,
%    containing the basis coefficients of the least-squares estimate
%    obtained by solving the equation
%     |A_11 ... A_r1|               |x_1|   |b_1|
%    (|...  ... ... | + lambda I) * |...| = |...|
%     |A_1r ... A_rr|               |x_r|   |b_r|
%    Each element of the block matrix A is a weighted sum of linear
%    mappings representing the projection/backprojection process.  lambda
%    is a regularization parameter; b is the set of backprojected image
%    sums at the appropriate sums, and with the appropriate weightings.
%
%    The equation is solved using the preconditioned conjugate gradient
%    method implemented by `conj_grad`.
%
%    mean_est_coeff: A vector of length `basis.count` containing the basis
%       coefficients of the least-squares estimate obtained by solving the
%       equation (A + lambda I)*x = b, where A is the linear mapping
%       represented by convolving with `mean_kernel_f`, lambda is the regular-
%       ization parameter `vols_wt_est_opt.regularizer`, and b is the sum of
%       backprojected images `vols_wt_b_coeff`, expressed in `basis`. The equa-
%       tion is solved using the preconditioned conjugate gradient method
%       implemented by `conj_grad`.
%    cg_info: A structure containing information about the conjugate gradient
%       method, such as residuals, objectives, etc. See the documentation of
%       `conj_grad` for more details.

% Author
%    Joakim Anden <janden@flatironinstitute.org>
%    Amit Halevi <ahalevi@princeton.edu>

function [mean_est_coeff, cg_info] = conj_grad_vols_wt(kermat_f, ...
    vols_wt_b_coeff, basis, precond_kermat_f, vols_wt_est_opt)

    if nargin < 4
        precond_kermat_f = [];
    end

    if nargin < 5 || isempty(vols_wt_est_opt)
        vols_wt_est_opt = struct();
    end

    kernel_fun = @(vols_coeff)( ...
        apply_vols_wt_kernel(vols_coeff, kermat_f, basis, vols_wt_est_opt));

    if ~isempty(precond_kermat_f)
        precond_fun = @(vols_coeff)( ...
            apply_vols_wt_kernel(vols_coeff, 1./precond_kermat_f, basis, ...
            vols_wt_est_opt));

        vols_wt_est_opt.preconditioner = precond_fun;
    end

    [mean_est_coeff, ~, cg_info] = ...
        conj_grad(kernel_fun, vols_wt_b_coeff(:), vols_wt_est_opt);

    mean_est_coeff = reshape(mean_est_coeff, [], size(vols_wt_b_coeff, 2));
end
