% IM_NORMALIZE_BG Normalize image background to zero mean and unit variance
%
% Usage
%    im = im_normalize_bg(im, bg_radius, do_ramp);
%
% Input
%    im: The image stack to process in the form of an L-by-L-by-n array.
%    bg_radius: The radius of the disk whose complement specifies the
%       background noise domain (default 1).
%    do_ramp: Fit a ramping background to the data and subtract. If false, a
%       constant background level is used instead (default true).
%
% Output
%    im: The same stack of images, but processed so their background noise has
%      zero mean and unit variance.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function im = im_normalize_bg(im, bg_radius, do_ramp)
    if nargin < 2 || isempty(bg_radius)
        bg_radius = 1;
    end

    if nargin < 3 || isempty(do_ramp)
        do_ramp = true;
    end

    [im, sz_roll] = unroll_dim(im, 3);

    g2d = grid_2d(size(im, 1));

    mask = find(g2d.r(:) > bg_radius);

    im = im_to_vec(im);

    if do_ramp
        design_mask = [g2d.x(mask) g2d.y(mask) ones(numel(mask), 1)];
        design_all = [g2d.x(:) g2d.y(:) ones(size(im, 1), 1)];

        coeff = design_mask\im(mask,:);

        im = im - design_all*coeff;
    end

    first_moment = mean(im(mask,:), 1);
    second_moment = mean(abs(im(mask,:)).^2, 1);

    bg_mean = first_moment;
    bg_std = sqrt(second_moment - first_moment.^2);

    im = bsxfun(@minus, im, bg_mean);
    im = bsxfun(@times, im, 1./bg_std);

    im = vec_to_im(im);

    im = roll_dim(im, sz_roll);
end
