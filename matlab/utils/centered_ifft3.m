% CENTERED_IFFT3 Calculate a centered, three-dimensional inverse FFT
%
% Usage
%    x = centered_ifft3(x_f);
%
% Input
%    x_f: The three-dimensional signal to be transformed. The inverse FFT is
%       only applied along the first three dimensions.
%
% Output
%    x: The centered inverse Fourier transform of x_f.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = centered_ifft3(x)
    x = ifftshift3(x);
    x = ifft3(x);
    x = fftshift3(x);
end
