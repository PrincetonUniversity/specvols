% SRC_MEAN_BACKWARD Apply adjoint mapping to source
%
% Usage
%    mean_b = src_mean_backward(src, mean_est_opt);
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
%    mean_b: The adjoint mapping applied to the images, averaged over the
%       whole dataset.
%
% See also
%    im_backward, src_mean_kernel

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function mean_b = src_mean_backward(src, mean_est_opt)
    if nargin < 2 || isempty(mean_est_opt)
        mean_est_opt = struct();
    end

    mean_est_opt = fill_struct(mean_est_opt, ...
        'precision', 'single', ...
        'batch_size', 512);

    params = src.params;

    L = src.L;
    n = src.n;

    if n ~= size(params.rots, 3)
        error('Number of images in source and parameters do not agree.');
    end

    batch_ct = ceil(n/mean_est_opt.batch_size);

    mean_b = zeros(L*ones(1, 3), mean_est_opt.precision);

    for batch = 1:batch_ct
        batch_s = (batch-1)*mean_est_opt.batch_size+1;
        batch_n = min(batch*mean_est_opt.batch_size, n)-batch_s+1;

        batch_mean_b = ...
            1/n*im_backward(src, src_image(src, batch_s, batch_n), batch_s);
        batch_mean_b = cast(batch_mean_b, mean_est_opt.precision);

        mean_b = mean_b + batch_mean_b;
    end
end
