% ROTATED_GRIDS Generate rotated Fourier grids in 3D from rotation matrices
%
% Usage
%    pts_rot = rotated_grids(N, rot_matrices);
%
% Input
%    N: The resolution of the desired grids.
%    rot_matrices: An array of size 3-by-3-by-K containing K rotation matrices.
%
% Output
%    pts_rot: A set of rotated Fourier grids in three dimensions as specified
%       by the rotation matrices. Frequencies are in the range [-pi, pi].

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function pts_rot = rotated_grids(L, rot_matrices)
    grid2d = grid_2d(L);

    num_rots = size(rot_matrices, 3);

    pts = pi*[grid2d.x(:) grid2d.y(:)];

    pts_rot = bsxfun(@times, rot_matrices(:,1,:), pts(:,1)') + ...
        bsxfun(@times, rot_matrices(:,2,:), pts(:,2)');

    pts_rot = reshape(pts_rot, [3 L*ones(1, 2) num_rots]);
end
