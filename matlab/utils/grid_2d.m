% GRID_2D Define an image grid
%
% Usage
%    grid = grid_2d(L);
%
% Input
%    L: The dimension of the desired grid.
%
% Output
%    mesh: A structure containing the fields:
%          - x, y: cartesian coordinates in an L-by-L grid, and
%          - r, phi: polar coordinates on the grid.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function grid = grid_2d(L)
    persistent p_grids;

    if size(p_grids, 1) >= L && ~isempty(p_grids{L})
        grid = p_grids{L};
        return;
    end

    grid1d = ceil([-L/2:L/2-1])/(L/2);

    % NOTE: Need to use ndgrid because meshgrid swaps x and y...
    [grid.x, grid.y] = ndgrid(grid1d, grid1d);

    [grid.phi, grid.r] = cart2pol(grid.x, grid.y);

    p_grids{L} = grid;
end
