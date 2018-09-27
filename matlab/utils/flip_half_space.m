% FLIP_HALF_SPACE Flip array to swap "positive" and "negative" half-spaces
%
% Usage
%    y = flip_half_space(x, dims);
%
% Input
%    x: The array to be flipped.
%    dims: A vector of dimensions along which to flip (default `1:ndims(x)`).
%
% Output
%    y: The array x, but flipped along the specified dimensions.
%
% See also
%    positive_half_space

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = flip_half_space(x, dims)
    if nargin < 2 || isempty(dims)
        dims = 1:numel(x);
    end

    idx_src.type = '()';
    idx_dst.type = '()';

    colon = {':'};

    idx_dst.subs = colon(ones(1, ndims(x)));
    idx_src.subs = colon(ones(1, ndims(x)));

    for dim = dims
        idx_dst.subs{dim} = 1:size(x, dim);
        idx_src.subs{dim} = [1 size(x, dim):-1:2];
    end

    x = mdim_ifftshift(x, dims);

    x = subsasgn(x, idx_dst, subsref(x, idx_src));

    x = mdim_fftshift(x, dims);
end
