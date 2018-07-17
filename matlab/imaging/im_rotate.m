% IM_ROTATE Rotate image
%
% Usage
%    im = im_rotate(im, angle);
%
% Input
%    im: An image of size L-by-L.
%    angle: An angle in the range 0 to 2*pi through which to rotate the image.
%
% Output
%    im: The original image, rotated clockwise by `angle`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function im = im_rotate(im, angle)
    L = size(im, 1);

    n = numel(angle);

    g2d = grid_2d(L);

    pts = pi*[g2d.x(:) g2d.y(:)]';

    R = angles_to_rots(cat(1, angle(:)', zeros(2, n)));
    R = R(1:2,1:2,:);

    pts_rot = zeros(2, L^2, n);
    for k = 1:n
        pts_rot(:,:,k) = R(:,:,k)*pts;
    end

    pts_rot = reshape(pts_rot, [2 L^2*n]);

    im_f = nufft2(im, pts_rot);

    im_f = reshape(im_f, [L*ones(1, 2) n]);

    if mod(L, 2) == 0
        im_f(1,:,:) = 0;
        im_f(:,1,:) = 0;
    end

    im = centered_ifft2(im_f);

    im = real(im);
end
