% VOL_FORWARD Apply forward imaging model to volume
%
% Usage
%    im = vol_forward(src, vol, s, n);
%
% Input
%    src: A source object from which to extract the imaging parameters. This
%       is typically obtained from `star_to_src` or `sim_to_src`.
%    vol: A volume of size L-by-L-by-L.
%    s: The first index of the parameters to use for the forward mapping.
%    n: The number of images to calculate.
%
% Output
%    im: The images obtained from `vol` by projecting, applying CTFs,
%       translating, and multiplying by the amplitude.
%
% See also
%    vol_project, im_backward

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function im = vol_forward(src, vol, s, n)
    params = src.params;

    L = size(vol, 1);

    idx = s:s+n-1;

    im = vol_project(vol, params.rots(:,:,idx));

    unique_filters = unique(params.filter_idx(idx));

    for k = unique_filters(:)'
        idx_k = find(params.filter_idx(idx) == k);

        im(:,:,idx_k) = im_filter(im(:,:,idx_k), params.filters(k));
    end

    im = im_translate(im, params.offsets(:,idx));

    im = bsxfun(@times, im, permute(params.amplitudes(idx), [1 3 2]));
end
