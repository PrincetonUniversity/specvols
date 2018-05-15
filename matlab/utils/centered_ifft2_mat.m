% CENTERED_IFFT2_MAT Calculate a centered, two-dimensional IFFT on a matrix
%
% Usage
%    X = centered_ifft2_mat(X_f);
%
% Input
%    X_f: The two-dimensional signal matrix to be transformed. The inverse FFT
%       is only applied along the first four dimensions.
%
% Output
%    X: The centered inverse Fourier transform of X along its two-dimensional
%       rows and columns.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function X = centered_ifft2_mat(X)
    dim_perm = [3 4 1 2 5:ndims(X)];

    X = centered_ifft2(X);
    X = conj(permute(centered_ifft2(permute(conj(X), dim_perm)), dim_perm));
end
