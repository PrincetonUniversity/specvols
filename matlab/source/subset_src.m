% SUBSET_SRC Extract subset of images from source
%
% Usage
%    src = subset_src(original_src, idx);
%
% Input
%    original_src: A source from which we want to define a subset source.
%    idx: A set of indices in `original_src` corresponding to the desired
%       subset.
%
% Output
%    src: A source containing the subset of images indicated by `idx`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function src = subset_src(original_src, idx)
    src = struct();

    src.type = src_type_subset();

    src.original_src = original_src;

    src.L = original_src.L;
    src.n = numel(idx);
    src.precision = original_src.precision;

    src.subset_idx = idx;

    src.params = struct();

    src.params.rots = original_src.params.rots(:,:,idx);
    src.params.filters = original_src.params.filters;
    src.params.filter_idx = original_src.params.filter_idx(:,idx);
    src.params.offsets = original_src.params.offsets(:,idx);
    src.params.amplitudes = original_src.params.amplitudes(:,idx);
end
