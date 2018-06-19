% ESTIMATE_NOISE Estimate noise variance from source
%
% Usage
%    noise_var_est = estimate_noise(src, noise_est_opt);
%
% Input
%    src: A source structure containing the images whose noise variance is to
%       be estimated.
%    noise_est_opt: A struct containing the fields:
%          - 'batch_size': The size of the batches in which to compute the
%             variance estimate (default 512).
%
% Output
%    noise_var_est: The estimated noise variance of the images.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function noise_var_est = estimate_noise(src, noise_est_opt)
    if nargin < 2 || isempty(noise_est_opt)
        noise_est_opt = struct();
    end

    noise_est_opt = fill_struct(noise_est_opt, ...
        'batch_size', 512);

    g2d = grid_2d(src.L);

    mask = find(g2d.r(:) >= 1);

    batch_ct = ceil(src.n/noise_est_opt.batch_size);

    first_moment = 0;
    second_moment = 0;

    for batch = 1:batch_ct
        batch_s = (batch-1)*noise_est_opt.batch_size+1;
        batch_n = min(batch*noise_est_opt.batch_size, src.n)-batch_s+1;

        im = src_image(src, batch_s, batch_n);

        im = im_to_vec(im);
        im_masked = im(mask,:);

        first_moment = first_moment + sum(im_masked(:))/(src.n*numel(mask));
        second_moment = second_moment + sum(abs(im_masked(:)).^2)/(src.n*numel(mask));
    end

    noise_var_est = second_moment-first_moment^2;
end
