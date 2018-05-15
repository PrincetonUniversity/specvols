% CENTERED_IFFT_MAT Calculate a centered, one-dimensional IFFT on a matrix
%
% Usage
%    X = centered_ifft_mat(X_f);
%
% Input
%    X_f: The one-dimensional signal matrix to be transformed. The inverse FFT
%       is only applied along the first and second dimensions.
%
% Output
%    X: The centered inverse Fourier transform of X along its rows and columns.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function X = centered_ifft_mat(X)
    dim_perm = [2 1 3:ndims(X)];

    X = centered_ifft(X);
    X = conj(permute(centered_ifft(permute(conj(X), dim_perm)), dim_perm));
end
