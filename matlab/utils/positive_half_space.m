% POSITIVE_HALF_SPACE Mask of "positive" half-space on a grid
%
% Usage
%    [mask, stat_mask] = positive_half_space(sz);
%
% Input
%    sz: Size of the grid.
%
% Output
%    mask: A mask of those coefficients corresponding to "positive" points on
%       the grid.
%    stat_mask: A mask of those coefficients which are invariant under
%       flipping from "positive" to "negative" half-space.
%
% See also
%    flip_half_space

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [mask, stat_mask] = positive_half_space(sz)
    for dim = 1:numel(sz)
        rngs{dim} = [0:floor(sz(dim)/2) -ceil(sz(dim)/2)+1:-1];
    end

    grids = cell(1, numel(sz));

    [grids{:}] = ndgrid(rngs{:});

    zero_mask = true(size(grids{1}));
    nonzero_mask = false(size(grids{1}));

    for dim = 1:numel(sz)
        nonzero_mask = nonzero_mask | ...
            (zero_mask&grids{dim}>0&grids{dim}<sz(dim)/2);
        zero_mask = zero_mask&(grids{dim}==0|grids{dim}==sz(dim)/2);
    end

    mask = zero_mask|nonzero_mask;
    stat_mask = zero_mask;

    mask = mdim_fftshift(mask, 1:ndims(sz));
    stat_mask = mdim_fftshift(stat_mask, 1:ndims(sz));
end
