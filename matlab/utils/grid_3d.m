% GRID_3D Define a volume grid
%
% Usage
%    grid = grid_3d(L);
%
% Input
%    L: The dimension of the desired grid.
%
% Output
%    grid: A structure containing the fields:
%          - x, y, z: cartesian coordinates in an L-by-L-by-L grid, and
%          - r, theta, phi: spherical coordinates on the grid.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function grid = grid_3d(L)
    persistent p_grids;

    if size(p_grids, 1) >= L && ~isempty(p_grids{L})
        grid = p_grides{L};
        return;
    end

    grid1d = ceil([-L/2:L/2-1])/(L/2);

    [grid.x, grid.y, grid.z] = ndgrid(grid1d, grid1d, grid1d);

    [grid.phi, grid.theta, grid.r] = cart2sph(grid.x, grid.y, grid.z);

    grid.theta = pi/2 - grid.theta;

    p_grids{L} = grid;
end
