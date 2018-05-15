% CENTERED_IFFT3_MAT Calculate a centered, three-dimensional IFFT on a matrix
%
% Usage
%    X = centered_ifft3_mat(X_f);
%
% Input
%    X_f: The three-dimensional signal matrix to be transformed. The inverse
%       FFT is only applied along the first six dimensions.
%
% Output
%    X: The centered inverse Fourier transform of X along its three-dimensional
%       rows and columns.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function X = centered_ifft3_mat(X)
    dim_perm = [4 5 6 1 2 3 7:ndims(X)];

    X = centered_ifft3(X);
    X = conj(permute(centered_ifft3(permute(conj(X), dim_perm)), dim_perm));
end
