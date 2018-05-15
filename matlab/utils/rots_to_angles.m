% ROTS_TO_ANGLES Calculate Euler angles from rotation matrices
%
% Usage
%    angles = rots_to_angles(rots);
%
% Input
%    rots: A 3-by-3-by-n array of rotation matrices concatenated
%       along the third dimension.
%
% Output
%    angles: The corresponding Euler angles.

% TODO: Specify which convention is used in the documentation.
% TODO: Speed up calculation by vectorizing atan2, acos calls.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function angles = rots_to_angles(rots)
    angles = zeros([3 size(rots, 3)]);

    for s = 1:size(rots, 3)
        [alpha, beta, gamma] = eulang(rots(:,:,s));
        angles(:,s) = [alpha beta gamma]';
    end
end

function [alpha, beta, gamma] = eulang(R)
    alpha = atan2(R(2,3), R(1,3));
    beta = acos(R(3,3));
    gamma = atan2(R(3,2), -R(3,1));
end
