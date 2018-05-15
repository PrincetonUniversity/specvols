% CENTERED_FFT3_MAT Calculate a centered, three-dimensional FFT on a matrix
%
% Usage
%    X_f = centered_fft3_mat(X);
%
% Input
%    X: The three-dimensional signal matrix to be transformed. The FFT is only
%       applied along the first six dimensions.
%
% Output
%    X_f: The centered Fourier transform of X along its three-dimensional rows
%       and columns.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function X = centered_fft3_mat(X)
    dim_perm = [4 5 6 1 2 3 7:ndims(X)];

    X = centered_fft3(X);
    X = conj(permute(centered_fft3(permute(conj(X), dim_perm)), dim_perm));
end
