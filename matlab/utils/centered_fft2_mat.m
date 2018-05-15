% CENTERED_FFT2_MAT Calculate a centered, two-dimensional FFT on a matrix
%
% Usage
%    X_f = centered_fft2_mat(X);
%
% Input
%    X: The two-dimensional signal matrix to be transformed. The FFT is only
%       applied along the first four dimensions.
%
% Output
%    X_f: The centered Fourier transform of X along its two-dimensional rows
%       and columns.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function X = centered_fft2_mat(X)
    dim_perm = [3 4 1 2 5:ndims(X)];

    X = centered_fft2(X);
    X = conj(permute(centered_fft2(permute(conj(X), dim_perm)), dim_perm));
end
