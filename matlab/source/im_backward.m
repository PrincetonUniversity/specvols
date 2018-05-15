% IM_BACKWARD Apply adjoint mapping to set of images
%
% Usage
%    vol = im_backward(src, im, s);
%
% Input
%    src: A source object from which to extract the imaging parameters. This
%       is typically obtained from `star_to_src` or `sim_to_src`.
%    im: An L-by-L-by-n array of images to which we wish to apply the adjoint
%       of the forward model.
%    s: The first index of the parameters to use for the adjoint mapping.
%
% Output
%    vol: An L-by-L-by-L volume containing the sum of the adjoint mappings
%       applied to the images.
%
% See also
%    im_backproject, vol_forward

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vol = im_backward(src, im, s)
    params = src.params;

    L = size(im, 1);
    n = size(im, 3);

    idx = s:s+n-1;

    im = bsxfun(@times, im, permute(params.amplitudes(idx), [1 3 2]));

    im = im_translate(im, -params.offsets(:,idx));

    unique_filters = unique(params.filter_idx(idx));

    for k = unique_filters(:)'
        idx_k = find(params.filter_idx(idx) == k);

        im(:,:,idx_k) = im_filter(im(:,:,idx_k), params.filters(k));
    end

    vol = im_backproject(im, params.rots(:,:,idx));
end
