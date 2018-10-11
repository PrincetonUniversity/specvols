% CACHE_SRC Cache whole source in memory
%
% Usage
%    src = cache_src(original_src);
%
% Input
%    original_src: The source to be cached.
%    cache_opt: A struct containing the fields:
%       - 'batch_size': The size of the batches used to cache the source
%          (default 512).
%
% Output
%    src: An in-memory array source containing the same images as
%       `original_src`.

function src = cache_src(original_src, cache_opt)
    if nargin < 2 || isempty(cache_opt)
        cache_opt = [];
    end

    cache_opt = fill_struct(cache_opt, ...
        'batch_size', 512);

    src = struct();

    src.type = src_type_array();

    src.L = original_src.L;
    src.n = original_src.n;
    src.precision = original_src.precision;

    src.images = zeros([src.L*ones(1, 2) src.n], src.precision);

    batch_ct = ceil(src.n/cache_opt.batch_size);

    for batch = 1:batch_ct
        batch_s = (batch-1)*cache_opt.batch_size+1;
        batch_n = min(batch*cache_opt.batch_size, src.n)-batch_s+1;

        batch_idx = batch_s:batch_s+batch_n-1;

        src.images(:,:,batch_idx) = src_image(original_src, batch_s, batch_n);
    end

    src.params = original_src.params;
end
