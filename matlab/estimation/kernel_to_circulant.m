% KERNEL_TO_CIRCULANT Convert 3D Fourier kernel into 3D circulant matrix
%
% Usage
%    A = kernel_to_circulant(kernel_f);
%
% Input
%    kernel_f: The Fourier transform of the convolution kernel in the form of
%       an L-by-L-by-L volume. The Fourier transform must be centered.
%
% Output
%    A: An six-dimensional circulant matrix of size L describing the convolu-
%       tion of a volume with the `kernel_f`.
%
% See also
%    kernel_to_toeplitz, conv_vol, src_mean_kernel, vec_to_vol, vol_to_vec

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function A = kernel_to_circulant(kernel_f)
    L = size(kernel_f, 1);

    A = kernel_to_toeplitz(kernel_f, L);
end
