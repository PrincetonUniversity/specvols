% ANGLES_TO_ROTS Calculate rotation matrices from Euler angles
%
% Usage
%    rots = angles_to_rots(angles);
%
% Input
%    angles: Euler angles in a 3-by-n matrix.
%
% Output
%    rots: Corresponding rotation matrices in a 3-by-3-by-n array.

% TODO: Specify which convention is used in the documentation.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function rots = angles_to_rots(angles)
    rots = zeros([3*ones(1, 2) size(angles, 2)]);

    for s = 1:size(angles, 2)
        rots(:,:,s) = erot(angles(:,s));
    end
end

function R = erot(angles)
    R = zrot(angles(1))*yrot(angles(2))*zrot(angles(3));
end

function R = zrot(theta)
    R = eye(3);
    R([1 2],[1 2]) = tworot(theta);
end

function R = yrot(theta)
    R = eye(3);
    R([3 1],[3 1]) = tworot(theta);
end

function R = tworot(theta)
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
end
