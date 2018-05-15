% CACHE_SRC Cache whole source in memory
%
% Usage
%    src = cache_src(original_src);
%
% Input
%    original_src: The source to be cached.
%
% Output
%    src: An in-memory array source containing the same images as
%       `original_src`.

function src = cache_src(original_src)
    src = struct();

    src.type = src_type_array();

    src.L = original_src.L;
    src.n = original_src.n;
    src.precision = original_src.precision;

    src.images = src_image(original_src, 1, src.n);

    src.params = original_src.params;
end
