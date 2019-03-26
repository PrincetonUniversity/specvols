% CENTERED_FFT3 Calculate a centered, three-dimensional FFT
%
% Usage
%    x_f = centered_fft3(x);
%
% Input
%    x: The three-dimensional signal to be transformed. The FFT is only applied
%       along the first three dimensions.
%
% Output
%    x_f: The centered Fourier transform of x.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = centered_fft3(x)
    x = ifftshift3(x);
    x = fft3(x);
    x = fftshift3(x);
end
