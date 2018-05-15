% APPLY_COVAR_KERNEL Apply the covariance kernel represented by convolution
%
% Usage
%    covar_coeff = apply_covar_kernel(covar_coeff, kernel_f, basis, ...
%       covar_est_opt);
%
% Input
%    covar_coeff: The matrix of covariance coefficients to be convolved,
%       represented in the basis as a `basis.count`-by-`basis.count` array.
%    kernel_f: The non-centered Fourier transform of the convolution kernel
%       representing the covariance projection-backprojection operator as
%       obtained from `src_covar_kernel`.
%    basis: A basis object corresponding to the basis used to store
%       `vol_coeff`.
%    covar_est_opt: An options structure. Currently no options are used.
%
% Output
%    covar_coeff: The result of evaluating `covar_coeff` in the given basis,
%       convolving with the kernel given by `kernel_f`, and backprojecting into
%       the basis.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function covar_coeff = apply_covar_kernel(covar_coeff, kernel_f, basis, ...
    covar_est_opt)

    covar = basis_mat_evaluate(basis, covar_coeff);

    covar = conv_volmat(covar, kernel_f);

    covar_coeff = basis_mat_evaluate_t(basis, covar);
end
