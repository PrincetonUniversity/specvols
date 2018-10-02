% SRC_VOLS_WT_KERMAT Calculate weighted mean estimation kernel matrix from
% source
%
% Usage
%    kermat_f = src_vols_wt_kermat(src, wts, mean_est_opt)
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    wts: a matrix of weights, r x n.
%    mean_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double' or
%             'single' (default).
%          - 'batch_size': The size of the batches in which to compute the
%             kernel (default 512).
%
% Output
%   kermat_f: A 2*L-by-2*L-by-2*L-by-r-by r array.  This is best considered
%       as an r x r matrix of volumes; each volume is a weighted mean
%       least-squares estimator kernel.  Convolution with each of these
%       kernels is equivalent to performing a projection/backprojection on
%       a volume, with the appropriate amplitude modifiers and CTF, and
%       also a weighting term; the r^2 volumes are each of pairwise
%       products between the weighting vectors given by the columns of wts.
%       Note that this is a non-centered Fourier transform, so the zero
%       frequency is found at index 1.
%
% See also
%    src_vols_wt_backward

% Authors
%    Joakim Anden <janden@flatironinstitute.org>
%    Amit Halevi <ahalevi@princeton.edu>

function kermat_f = src_vols_wt_kermat_bsx(src, wts, mean_est_opt)
    if nargin < 2 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    params = src.params;

    L = src.L;
    n = size(params.rots, 3);
    r = size(wts, 1);

    kermat = zeros([2*L 2*L 2*L r r], mean_est_opt.precision);

    mean_est_opt = fill_struct(mean_est_opt, ...
        'precision', 'single', ...
        'batch_size', 512);

    batch_ct = ceil(n/mean_est_opt.batch_size);

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

    kermat_f = zeros(2*L*ones(1, 3), mean_est_opt.precision);

    %final clean-up
    % Ensure symmetric kernel.
    for j = 1:r
        for k = 1:r
            kermat(1,:,:,j,k) = 0;
            kermat(:,1,:,j,k) = 0;
            kermat(:,:,1,j,k) = 0;

            % Compute non-centered Fourier transform.
            kermat(:,:,:,j,k) = mdim_ifftshift(kermat(:,:,:,j,k), 1:3);
            kermat_f(:,:,:,j,k) = fft3(kermat(:,:,:,j,k));
        end
    end
    % Kernel is always symmetric in spatial domain and therefore real in
    % Fourier.
    kermat_f = real(kermat_f);
end
