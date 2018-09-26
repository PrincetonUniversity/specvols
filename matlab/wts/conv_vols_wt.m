% CONV_VOL Convolve volume(s) with kernel(s)
%
% Usage
%    vol = conv_vol(vol, kermat_f);
%
% Input
%    vol: An N-by-N-by-N-by-... array of volumes to be convolved.
%    kermat_f: The Fourier transform of the cubic convolution kernel. Must be
%       larger than vol in the first three dimensions. This is a non-centered
%       Fourier transform. That is, the zero frequency is found at index 1.
%
% Output
%    vol: The original volumes convolved by the kernel with the same dimensions
%       as before.

% Author
%    Joakim Anden <janden@flatironinstitute.org>
%    Amit Halevi <ahalvei@princeton.edu>

function x = conv_vols_wt(x, kermat_f)
    N = size(x, 1);

    sz_x = size(x);
    sz_kermat = size(kermat_f);

    if any(2*sz_x(1:3) ~= sz_kermat(1:3)) && any(sz_x(1:3) ~= sz_kermat(1:3))
        error('Inconsistent spatial dimensions.');
    end
% 
    if ~isempty(sz_x(4:end)) && ~isempty(sz_kermat(4:end)) && ...
        (numel(sz_x) ~= numel(sz_kermat) || any(sz_x(4:end) ~= sz_kermat(4:end)))
        error('Inconsistent indexing dimensions.');
        % TODO: Talk about singletons dimensions.
    end

    if any(sz_x(1:3) ~= N)
        error('Volumes in `x` must be cubic.');
    end

    [x, sz_roll] = unroll_dim(x, 4);
    [kermat_f, ~] = unroll_dim(kermat_f, 4);

    N_ker = size(kermat_f, 1);

    x = fft(x, N_ker, 1);
    x = fft(x, N_ker, 2);
    x = fft(x, N_ker, 3);

    x = bsxfun(@times, x, kermat_f);

    x = ifft(x, [], 1);
    x = x(1:N,:,:,:);
    x = ifft(x, [], 2);
    x = x(:,1:N,:,:);
    x = ifft(x, [], 3);
    x = x(:,:,1:N,:);

    x = real(x);

    x = roll_dim(x, sz_roll);
end
