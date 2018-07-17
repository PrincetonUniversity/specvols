% VOL_ROTATE Rotate volume
%
% Usage
%    vol = vol_rotate(vol, rot_matrix);
%
% Input
%    vol: A volume of size L-by-L-by-L.
%    rot_matrix: A rotation matrix of size 3-by-3 through which to rotate the
%       volume.
%
% Output
%    vol: The original volume, rotated by `rot_matrix`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vol = vol_rotate(vol, rot_matrix)
    L = size(vol, 1);

    n = size(rot_matrix, 3);

    g3d = grid_3d(L);

    pts = pi*[g3d.x(:) g3d.y(:) g3d.z(:)]';

    pts_rot = zeros(3, L^3, n);
    for k = 1:n
        pts_rot(:,:,k) = rot_matrix(:,:,k)*pts;
    end

    pts_rot = reshape(pts_rot, [3 L^3*n]);

    vol_f = nufft3(vol, pts_rot);

    vol_f = reshape(vol_f, [L*ones(1, 3) n]);

    if mod(L, 2) == 0
        vol_f(1,:,:,:) = 0;
        vol_f(:,1,:,:) = 0;
        vol_f(:,:,1,:) = 0;
    end

    vol = centered_ifft3(vol_f);

    vol = real(vol);
end
