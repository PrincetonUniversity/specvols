% CONV_VOL_BKWD Convolve volume with kernels
%
% Usage
%    vol = conv_vol(vol, kernels_f);
%
% Input
%    vol: An N-by-N-by-N-by-r array of volumes to be convolved.
%    vol: A volume to be convolved by several kernels
%    kernels_f: a 2N-by-2N-by-2N-by... array of kernels with which to
%    convolve the volume col.  Note that these kernels are all non-centered
%    FTs - i.e., zero frequency is at index 1.
%
% Output
%    vols: The original volume convolved by each of the kernels, an array
%    of size N-by-N-by-N-by-r

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = conv_vol_bkwd(x, kernels_f)
    N = size(x, 1);

    [x, sz_roll] = unroll_dim(x, 4);

%     if any([size(x, 1) size(x, 2) size(x, 3)] ~= N)
%         error('Volumes in `x` must be cubic.');
%     end

%     is_singleton = (numel(size(x)) == 3);
    is_singleton = 0;
    N_ker = size(kernels_f, 1);

%     if any(size(kernel_f) ~= N_ker)
%         error('Convolution kernel `kernel_f` must be cubic.');
%     end

    if is_singleton
        x = fftn(x, N_ker*ones(1, 3));
    else
        x = fft(x, N_ker, 1);
        x = fft(x, N_ker, 2);
        x = fft(x, N_ker, 3);
    end

%     x = bsxfun(@times, x, kernel_f);

    x = bsxfun(@times, kernels_f, x);

    if is_singleton
        if ~isoctave
            x = ifftn(x, [], 'symmetric');
        else
            x = ifftn(x);
        end

        x = x(1:N,1:N,1:N,:);
    else
        x = ifft(x, [], 1);
        x = x(1:N,:,:,:);
        x = ifft(x, [], 2);
        x = x(:,1:N,:,:);
        x = ifft(x, [], 3);
        x = x(:,:,1:N,:);
    end

    x = real(x);

    x = roll_dim(x, sz_roll);
end
