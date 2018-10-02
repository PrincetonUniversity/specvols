% SRC_MEAN_KERNEL Calculate mean estimation kernel from source
%
% Usage
%    mean_kernel_f = src_mean_kernel(src, mean_est_opt);
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    mean_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double' or
%             'single' (default).
%          - 'batch_size': The size of the batches in which to compute the
%             kernel (default 512).
%
% Output
%    mean_kernel_f: A 2*L-by-2*L-by-2*L array containing the non-centered
%       Fourier transform of the mean least-squares estimator kernel. Convol-
%       ving a volume with this kernel is equal to projecting and backproject-
%       ing that volume in each of the projection directions (with the appro-
%       priate amplitude multipliers and CTFs) and averaging over the whole
%       dataset. Note that this is a non-centered Fourier transform, so the
%       zero frequency is found at index 1.
%
% See also
%    src_mean_backward

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function mean_kernel_f = src_mean_kernel(src, mean_est_opt)
    if nargin < 2 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    params = src.params;

    L = src.L;
    n = size(params.rots, 3);

    mean_est_opt = fill_struct(mean_est_opt, ...
        'precision', 'single', ...
        'batch_size', 512);

    batch_ct = ceil(n/mean_est_opt.batch_size);

    mean_kernel = zeros(2*L*ones(1, 3), mean_est_opt.precision);

    filters_f = eval_filter_grid(params.filters, L);
    sq_filters_f = abs(filters_f).^2;

    sq_filters_f = cast(sq_filters_f, mean_est_opt.precision);

    for batch = 1:batch_ct
        batch_s = (batch-1)*mean_est_opt.batch_size+1;
        batch_n = min(batch*mean_est_opt.batch_size, n)-batch_s+1;

        batch_idx = batch_s:batch_s+batch_n-1;

        pts_rot = rotated_grids(L, params.rots(:,:,batch_idx));

        weights = sq_filters_f(:,:,params.filter_idx(batch_idx));
        weights = bsxfun(@times, weights, ...
            reshape(params.amplitudes(batch_idx).^2, [1 1 batch_n]));

        if mod(L, 2) == 0
            weights(1,:,:) = 0;
            weights(:,1,:) = 0;
        end

        weights = im_to_vec(weights);

        pts_rot = reshape(pts_rot, 3, L^2*batch_n);
        weights = reshape(weights, L^2*batch_n, 1);

        mean_kernel = mean_kernel + ...
            1/(n*L^4)*real(anufft3(weights, pts_rot, 2*L*ones(1, 3)));
    end

    % Ensure symmetric kernel.
    mean_kernel(1,:,:) = 0;
    mean_kernel(:,1,:) = 0;
    mean_kernel(:,:,1) = 0;

    % Compute non-centered Fourier transform.
    mean_kernel = mdim_ifftshift(mean_kernel, 1:3);
    mean_kernel_f = fft3(mean_kernel);

    % Kernel is always symmetric in spatial domain and therefore real in
    % Fourier.
    mean_kernel_f = real(mean_kernel_f);
end
