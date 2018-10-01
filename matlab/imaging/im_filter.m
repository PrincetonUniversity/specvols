% IM_FILTER Apply filter to an image
%
% Usage
%    im = im_filter(x, filter);
%
% Input
%    x: Original array of size L-by-L-by-... whose first two dimensions
%       consist of images to be filtered.
%    filter: A filter object (see `eval_filter` for more information).
%
% Output
%    x: The array of filtered images.
%
% Note
%    The filtering is performed in a periodic manner, so for point spread
%    functions of large support, there may be boundary artifacts.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function im = im_filter(im, filter)
    if filter.type == filter_type_scalar()
        if filter.value ~= 1
            im = filter.value*im;
        end

        return;
    end

    L = size(im, 1);

    [im, sz_roll] = unroll_dim(im, 3);

    filter_vals = eval_filter_grid(filter, L);

    im_f = centered_fft2(im);

    im_f = bsxfun(@times, im_f, filter_vals);

    im = centered_ifft2(im_f);

    im = real(im);

    im = roll_dim(im, sz_roll);
end
