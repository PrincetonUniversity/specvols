% IM_BACKWARD_WT Apply adjoint mapping to set of images with weights
%
% Usage
%    vol = im_backward_wt(src, im, wts, s)
%
% Input
%    src: A source object from which to extract the imaging parameters. This
%       is typically obtained from `star_to_src` or `sim_to_src`.
%    im: An L-by-L-by-n array of images to which we wish to apply the adjoint
%       of the forward model.
%    wts: the set of weights to apply.  This is the entire set of weights
%    for all images, not just those in the current batch - i.e., it is a
%    vector of length equal to the entire number of images being used for
%    reconstruction, not equal to n (the batch size)
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
%    Amit Halevi <ahalevi@princeton.edu>

function vol = im_backward_wt_multi(src, im, wts, s)
    params = src.params;

    L = size(im, 1);
    n = size(im, 3);

    idx = s:s+n-1;

    im = bsxfun(@times, im, permute(params.amplitudes(idx), [1 3 2]));

    im = im_translate(im, -params.offsets(:,idx));

    unique_filters = unique(params.filter_idx(idx));

    for k = unique_filters(:)'
        idx_k = find(params.filter_idx(idx) == k);

        %There might be a "better" way of doing this - say by multiplying
        %the filters, but this works well enough for now
%         im(:,:,idx_k) = bsxfun(@times,im_filter(im(:,:,idx_k), ...
%             params.filters(k)),reshape(wts(idx_k),[1 1 numel(idx_k)]));
        im(:,:,idx_k) = im_filter(im(:,:,idx_k), params.filters(k));
    end

    vol = im_backproject(im, params.rots(:,:,idx),wts);
end
