% DOWNSAMPLE_SRC Blur and downsample source
%
% Usage
%    src = downsample_src(original_src, L);
%
% Input
%    original_src: A source whose images are to be downsampled.
%    L: The desired resolution of those images. Must be smaller than `src.L`.
%
% Output
%    src: A source containing the images as the source, but downsampled to
%       resolution `L`.
%
% See also
%    im_downsample

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function src = downsample_src(original_src, L)
    src = struct();

    src.type = src_type_downsampled();

    src.original_src = original_src;

    src.L = L;
    src.n = original_src.n;
    src.precision = original_src.precision;

    ds_factor = original_src.L/src.L;

    src.params = original_src.params;

    src.params.offsets = src.params.offsets/ds_factor;

    for k = 1:numel(src.params.filters)
        src.params.filters(k) = scale_filter(src.params.filters(k), ds_factor);
    end
end

function filter = scale_filter(filter, c)
    if filter.type == filter_type_scalar()
        return;
    elseif filter.type == filter_type_ctf()
        filter.ctf_params.pixel_size = filter.ctf_params.pixel_size*c;
        return;
    end
end
