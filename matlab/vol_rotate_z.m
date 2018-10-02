% VOL_ROTATE_Z Rotate volume around z-axis
%
% Usage
%    vol_rot = vol_rotate_z(vol, theta);
%
% Input
%    vol: An N-by-N-by-N volume array to be rotated.
%    theta: A vector of M angles with which to rotate the volume (in radians).
%
% Output
%    vol_rot: An N-by-N-by-N-by-M array consisting of the volume rotated through
%       the various angles.

function vol_rot = vol_rotate_z(vol, theta)
    vol_rot = zeros([size(vol) numel(theta)], class(vol));

    theta_deg = theta/pi*180;

    rotate_image = @(im, angle)(imrotate(im, angle, 'bilinear', 'crop'));

    for k = 1:numel(theta)
        vol_rot(:,:,:,k) = submatfun(@(im)(rotate_image(im, theta_deg(k))), vol, 3);
    end
end
