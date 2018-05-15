% CONV_VOLMAT Convolve volume matrix with kernel
%
% Usage
%    X = conv_volmat(X, kernel_f);
%
% Input
%    X: An N-by-...-by-N (6 dimensions) volume matrix to be convolved.
%    kernel_f: The Fourier transform of the cubic matrix convolution kernel.
%        Must be larger than volmat. This must be a non-centered Fourier
%        transform.
%
% Output
%    X: The original volume matrix convolved by the kernel with the same
%       dimensions as before.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function X = conv_volmat(X, kernel_f)
    sz_volmat = size(X);
    N = sz_volmat(1);

    if any(sz_volmat(1:6) ~= N)
        error('Volume matrix must be cubic and square.');
    end

    % TODO: Deal with rolled dimensions

    is_singleton = (numel(sz_volmat) == 6);

    N_ker = size(kernel_f, 1);

    if any(size(kernel_f) ~= N_ker)
        error('Convolution kernel must be cubic and square.');
    end

    if ~isoctave()
        % NOTE: Order is important here. It's about 20% faster to run from
        % 1 through 6 compared with 6 through 1.
        X = fft(X, N_ker, 1);
        X = fft(X, N_ker, 2);
        X = fft(X, N_ker, 3);
        X = fft(X, N_ker, 4);
        X = fft(X, N_ker, 5);
        X = fft(X, N_ker, 6);

        X = bsxfun(@times, X, kernel_f);

        % NOTE: Again, the order here is important. Also, we can't use
        % 'symmetric' option here since the partial IFFTs are not real.
        X = ifft(X, [], 6);
        X = X(:,:,:,:,:,1:N,:);
        X = ifft(X, [], 5);
        X = X(:,:,:,:,1:N,:,:);
        X = ifft(X, [], 4);
        X = X(:,:,:,1:N,:,:,:);
        X = ifft(X, [], 3);
        X = X(:,:,1:N,:,:,:,:);
        X = ifft(X, [], 2);
        X = X(:,1:N,:,:,:,:,:);
        X = ifft(X, [], 1);
        X = X(1:N,:,:,:,:,:,:);
    else
        if is_singleton
            X = fftn(X, N_ker*ones(1, 6));
        else
            X = fft(X, N_ker, 1);
            X = fft(X, N_ker, 2);
            X = fft(X, N_ker, 3);
            X = fft(X, N_ker, 4);
            X = fft(X, N_ker, 5);
            X = fft(X, N_ker, 6);
        end

        X = bsxfun(@times, X, kernel_f);

        if is_singleton
            X = ifftn(X);
            X = X(1:N,1:N,1:N,1:N,1:N,1:N);
        else
            X = ifft(X, [], 6);
            X = X(:,:,:,:,:,1:N,:);
            X = ifft(X, [], 5);
            X = X(:,:,:,:,1:N,:,:);
            X = ifft(X, [], 4);
            X = X(:,:,:,1:N,:,:,:);
            X = ifft(X, [], 3);
            X = X(:,:,1:N,:,:,:,:);
            X = ifft(X, [], 2);
            X = X(:,1:N,:,:,:,:,:);
            X = ifft(X, [], 1);
            X = X(1:N,:,:,:,:,:,:);
        end
    end

    X = real(X);
end

