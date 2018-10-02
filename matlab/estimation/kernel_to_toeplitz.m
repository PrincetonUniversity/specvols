% KERNEL_TO_TOEPLITZ Convert 3D Fourier kernel into 3D Toeplitz matrix
%
% Usage
%    A = kernel_to_toeplitz(kernel_f);
%
% Input
%    kernel_f: The Fourier transform of the convolution kernel in the form of
%       2*L-by-2*L-by-2*L volume. The Fourier transform must be centered.
%
% Output
%    A: An six-dimensional Toeplitz matrix of size L describing the convolu-
%       tion of a volume with the `kernel_f`.
%
% See also
%    conv_vol, src_mean_kernel, vec_to_vol, vol_to_vec

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function A = kernel_to_toeplitz(kernel_f)
    L = size(kernel_f, 1)/2;

    A = eye(L^3, class(kernel_f));

    for k = 1:L^3
        A(:,k) = vol_to_vec(conv_vol(vec_to_vol(A(:,k)), kernel_f));
    end

    A = vecmat_to_volmat(A);
end
