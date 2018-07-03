% SRC_VOLS_WT_BACKWARD Apply weighted adjoint mapping to source with
% weights
%
% Usage
%    vols_wt_b_coeff = src_vols_wt_backward(src, basis, wts, mean_est_opt);
%
% Input
%    src: A source structure containing the images and imaging parameters.
%       This is typically obtained from `star_to_src` or `sim_to_src`.
%    basis: A basis object used for representing the volumes.
%    wts: a vector of weights to apply to the respective images when
%    finding the mean of backprojected images
%    mean_est_opt: A struct containing the fields:
%          - 'precision': The precision of the kernel. Either 'double' or
%             'single' (default).
%          - 'batch_size': The size of the batches in which to compute the
%             kernel (default 512).
%
% Output
%    vols_wt_b_coeff: A set of volumes, each of which is a copy of the adjoint
%       mapping applied to the images, averaged with weights given by the
%       respective columns of wts, over the whole dataset and expressed as
%       coefficients of "basis".  A basis.count x r matrix
%
% See also
%    im_backward_wt, src_vols_wt_kermat

% Author
%    Joakim Anden <janden@flatironinstitute.org>
%    Amit Halevi <ahalevi@princeton.edu>

function vols_wt_b_coeff = src_vols_wt_backward(src, basis, wts, mean_est_opt)
    if nargin < 3 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    mean_est_opt = fill_struct(mean_est_opt, ...
        'precision', 'single', ...
        'batch_size', 512);

    params = src.params;

    L = src.L;
    n = src.n;
    r = size(wts,1);

    if n ~= size(params.rots, 3)
        error('Number of images in source and parameters do not agree.');
    end

    batch_ct = ceil(n/mean_est_opt.batch_size);

    vols_wt_b = zeros([L L L r], mean_est_opt.precision);

    for k = 1:r
        for batch = 1:batch_ct
            batch_s = (batch-1)*mean_est_opt.batch_size+1;
            batch_n = min(batch*mean_est_opt.batch_size, n)-batch_s+1;

            batch_vols_wt_b = 1/n * im_backward_wt(src, src_image(src, ...
                batch_s, batch_n), wts(k,:), batch_s);
            batch_vols_wt_b = cast(batch_vols_wt_b, mean_est_opt.precision);

            vols_wt_b(:,:,:,k) = vols_wt_b(:,:,:,k) + batch_vols_wt_b;
        end
    end

    vols_wt_b_coeff = basis_evaluate_t(basis, vols_wt_b);
end
