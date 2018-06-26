% ESTIMATE_NOISE Estimate noise power spectrum from source
%
% Usage
%    noise_psd_est = estimate_noise(src, noise_est_opt);
%
% Input
%    src: A source structure containing the images whose noise power spectrum
%       is to be estimated.
%    noise_est_opt: A struct containing the fields:
%          - 'batch_size': The size of the batches in which to compute the
%             variance estimate (default 512).
%          - 'noise_type': The type of noise to estimate. This is either
%             'white', in which case only the variance is estimated, or
%             'anisotropic', in which an anisotropic power spectral
%             distribution is estimated (default 'white').
%
% Output
%    noise_psd_est: The estimated noise power spectral distribution (PSD) of
%       the images in the form of a filter object.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function noise_psd_est = estimate_noise(src, noise_est_opt)
    if nargin < 2 || isempty(noise_est_opt)
        noise_est_opt = struct();
    end

    noise_est_opt = fill_struct(noise_est_opt, ...
        'batch_size', 512, ...
        'noise_type', 'white');

    if strcmp(noise_est_opt.noise_type, 'white')
        noise_var_est = estimate_noise_white(src, noise_est_opt);
        noise_psd_est = scalar_filter(noise_var_est, 2);
    elseif strcmp(noise_est_opt.noise_type, 'anisotropic')
        noise_psd_est = estimate_noise_psd(src, noise_est_opt);
        noise_psd_est = array_filter(noise_psd_est);
    else
        error('Invalid `noise_type`.');
    end
end

function noise_psd_est = estimate_noise_psd(src, noise_est_opt)
    g2d = grid_2d(src.L);

    mask = (g2d.r(:) >= 1);

    batch_ct = ceil(src.n/noise_est_opt.batch_size);

    mean_est = 0;

    noise_psd_est = zeros(src.L*ones(1, 2), src.precision);

    for batch = 1:batch_ct
        batch_s = (batch-1)*noise_est_opt.batch_size+1;
        batch_n = min(batch*noise_est_opt.batch_size, src.n)-batch_s+1;

        im = src_image(src, batch_s, batch_n);

        im_masked = im_to_vec(im);
        im_masked(~mask,:) = 0;
        im_masked = vec_to_im(im_masked);

        mean_est = mean_est + sum(im_masked(:))/(src.n*sum(mask));

        im_masked_f = centered_fft2(im_masked);

        noise_psd_est = noise_psd_est + sum(abs(im_masked_f).^2, 3)/(src.n*sum(mask));
    end

    mid = floor(src.L/2)+1;
    noise_psd_est(mid,mid) = noise_psd_est(mid,mid) - mean_est^2;
end

function noise_var_est = estimate_noise_white(src, noise_est_opt)
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
