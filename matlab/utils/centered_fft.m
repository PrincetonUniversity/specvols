% CENTERED_FFT Calculate a centered, one-dimensional FFT
%
% Usage
%    x_f = centered_fft(x);
%
% Input
%    x: The one-dimensional signal to be transformed. The FFT is only applied
%       along the first dimension.
%
% Output
%    x_f: The centered Fourier transform of x.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = centered_fft(x)
    x = ifftshift(x, 1);
    x = fft(x, [], 1);
    x = fftshift(x, 1);
end
