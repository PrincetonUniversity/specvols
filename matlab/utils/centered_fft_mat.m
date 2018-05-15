% CENTERED_FFT_MAT Calculate a centered, one-dimensional FFT on a matrix
%
% Usage
%    X_f = centered_fft_mat(X);
%
% Input
%    X: The one-dimensional signal matrix to be transformed. The FFT is only
%       applied along the first and second dimensions.
%
% Output
%    X_f: The centered Fourier transform of X along its rows and columns.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function X = centered_fft_mat(X)
    dim_perm = [2 1 3:ndims(X)];

    X = centered_fft(X);
    X = conj(permute(centered_fft(permute(conj(X), dim_perm)), dim_perm));
end
