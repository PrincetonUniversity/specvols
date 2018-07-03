% CONV_VOL Convolve volume with kernel
%
% Usage
%    vol = conv_vol(vol, kernel_f);
%
% Input
%    vol: An N-by-N-by-N-by-... array of volumes to be convolved.
%    kernel_f: The Fourier transform of the cubic convolution kernel. Must be
%       larger than vol in the first three dimensions. This is a non-centered
%       Fourier transform. That is, the zero frequency is found at index 1.
%
% Output
%    vol: The original volumes convolved by the kernel with the same dimensions
%       as before.

% Author
%    Joakim Anden <janden@flatironinstitute.org>
%    Amit Halevi <ahalvei@princeton.edu>

function x = conv_vol_wt(x, kernel_f)
    N = size(x, 1);

    sz_x = size(x);
    sz_kernel = size(kernel_f);

    if any(2*sz_x(1:3) ~= sz_kernel(1:3))
        error('Inconsistent spatial dimensions.');
    end

    if ~isempty(sz_x(4:end)) && ~isempty(sz_kernel(4:end)) && ...
        (numel(sz_x) ~= numel(sz_kernel) || any(sz_x(4:end) ~= sz_kernel(4:end)))
        error('Inconsistent indexing dimensions.');
        % TODO: Talk about singletons dimensions.
    end

    if any(sz_x(1:3) ~= N)
        error('Volumes in `x` must be cubic.');
    end

    [x, sz_roll] = unroll_dim(x, 4);
    [kernel_f, ~] = unroll_dim(kernel_f, 4);

    N_ker = size(kernel_f, 1);

    x = fft(x, N_ker, 1);
    x = fft(x, N_ker, 2);
    x = fft(x, N_ker, 3);

    x = bsxfun(@times, x, kernel_f);

    x = ifft(x, [], 1);
    x = x(1:N,:,:,:);
    x = ifft(x, [], 2);
    x = x(:,1:N,:,:);
    x = ifft(x, [], 3);
    x = x(:,:,1:N,:);

    x = real(x);

    x = roll_dim(x, sz_roll);
end
