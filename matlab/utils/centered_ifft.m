% CENTERED_IFFT Calculate a centered, one-dimensional inverse FFT
%
% Usage
%    x = centered_ifft(x_f);
%
% Input
%    x_f: The one-dimensional signal to be transformed. The inverse FFT is
%       only applied along the first dimension.
%
% Output
%    x: The centered inverse Fourier transform of x_f.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = centered_ifft(x)
    x = ifftshift(x, 1);
    x = ifft(x, [], 1);
    x = fftshift(x, 1);
end
