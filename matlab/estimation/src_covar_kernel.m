% SRC_COVAR_KERNEL Calculate covariance estimation kernel from source
%
% Usage
%    covar_kernel_f = src_covar_kernel(src, covar_est_opt);
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    covar_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double' or
%             'single' (default).
%          - 'batch_size': The size of the batches in which to compute the
%             kernel (default 512).
%
% Output
%    covar_kernel_f: A 2*L-by-2*L-by-2*L-by-2*L-by-2*L-by-2*L array containing
%       the non-centered Fourier transform of the covariance least-squares
%       estimator kernel. Convolving a volume matrix with this kernel is equal
%       to projecting and backprojecting that matrix along both rows and
%       columns in each projection direction (with appropriate amplitude
%       multipliers, CTFs, and shifts), and averaging over the whole dataset.
%       Note that this is a non-centered Fourier transform, so the zero fre-
%       quency is found at index 1.
%
% See also
%    src_covar_backward

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function covar_kernel_f = src_covar_kernel(src, covar_est_opt)
    if nargin < 2 || isempty(covar_est_opt)
        covar_est_opt = struct();
    end

    covar_est_opt = fill_struct(covar_est_opt, ...
        'precision', 'single', ...
        'batch_size', 512);

    covar_kernel = zeros(2*src.L*ones(1, 6), covar_est_opt.precision);

    filters_f = eval_filter_grid(src.params.filters, src.L);
    sq_filters_f = abs(filters_f).^2;

    sq_filters_f = cast(sq_filters_f, covar_est_opt.precision);

    for batch = 1:ceil(src.n/covar_est_opt.batch_size)
        batch_s = (batch-1)*covar_est_opt.batch_size+1;
        batch_n = min(batch*covar_est_opt.batch_size, src.n)-batch_s+1;

        batch_idx = batch_s:batch_s+batch_n-1;

        pts_rot = rotated_grids(src.L, src.params.rots(:,:,batch_idx));

        weights = sq_filters_f(:,:,src.params.filter_idx(batch_idx));
        weights = bsxfun(@times, weights, ...
            reshape(src.params.amplitudes(batch_idx).^2, [1 1 batch_n]));

        if mod(src.L, 2) == 0
            weights(1,:,:) = 0;
            weights(:,1,:) = 0;
        end

        weights = im_to_vec(weights);

        pts_rot = reshape(pts_rot, [3, src.L^2, batch_n]);
        weights = reshape(weights, [src.L^2, batch_n]);

        factors = zeros([2*src.L*ones(1, 3) batch_n], ...
            covar_est_opt.precision);

        for s = 1:batch_n
            factors(:,:,:,s) = real(anufft3(weights(:,s), ...
                pts_rot(:,:,s), 2*src.L*ones(1, 3)));
        end

        factors = vol_to_vec(factors);

        covar_kernel = covar_kernel + ...
            1/(src.n*src.L^8)*vecmat_to_volmat(factors*factors');
    end

    % Ensure symmetric kernel.
    covar_kernel(1,:,:,:,:,:) = 0;
    covar_kernel(:,1,:,:,:,:) = 0;
    covar_kernel(:,:,1,:,:,:) = 0;
    covar_kernel(:,:,:,1,:,:) = 0;
    covar_kernel(:,:,:,:,1,:) = 0;
    covar_kernel(:,:,:,:,:,1) = 0;

    % Compute non-centered Fourier transform.
    covar_kernel = mdim_ifftshift(covar_kernel, 1:6);
    covar_kernel_f = fftn(covar_kernel);

    % Kernel is always symmetric in spatial domain and therefore real in
    % Fourier.
    covar_kernel_f = real(covar_kernel_f);
end
