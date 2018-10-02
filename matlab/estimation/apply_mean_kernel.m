% APPLY_MEAN_KERNEL Applies the mean kernel represented by convolution
%
% Usage
%    vol_coeff = apply_mean_kernel(vol_coeff, kernel_f, basis, mean_est_opt);
%
% Input
%    vol_coeff: The volume to be convolved, stored in the basis coefficients.
%    kernel_f: The non-centered Fourier transform of the convolution kernel
%       representing the mean projection-backprojection operator as obtained
%       from `src_mean_kernel`.
%    basis: A basis object corresponding to the basis used to store
%       `vol_coeff`.
%    mean_est_opt: An options structure. Currently no options are used.
%
% Output
%    vol_coeff: The result of evaluating `vol_coeff` in the given basis,
%       convolving with the kernel given by `kernel_f`, and backprojecting into
%       the basis.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vol_coeff = apply_mean_kernel(vol_coeff, kernel_f, basis, ...
    mean_est_opt)

    vol = basis_evaluate(basis, vol_coeff);

    vol = conv_vol_bkwd(vol, kernel_f);

    vol_coeff = basis_evaluate_t(basis, vol);
end
