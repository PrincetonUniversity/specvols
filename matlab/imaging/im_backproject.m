% IM_BACKPROJECT Backproject images along rotation
%
% Usage
%    vol = im_backproject(im, rot_matrices);
%
% Input
%    im: An L-by-L-by-n array of images to backproject.
%    rot_matrices: An 3-by-3-by-n array of rotation matrices corresponding to
%       viewing directions.
%
% Output
%    vol: An L-by-L-by-L volumes corresponding to the sum of the backprojected
%       images.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vol = im_backproject(im, rot_matrices)
    L = size(im, 1);

    if size(im, 2) ~= L
        error('im must be of the form L-by-L-by-K');
    end

    n = size(im, 3);

    if size(rot_matrices, 3) ~= n
        error(['The number of rotation matrices must match the number ' ...
            'of images.']);
    end

    [real_mask, stat_mask] = positive_half_space(L*ones(1, 2));

    pts_rot = rotated_grids(L, rot_matrices, real_mask);

    pts_rot = reshape(pts_rot, [3 sum(real_mask(:))*n]);

    im_f = 1/L^2*centered_fft2(im);

    if mod(L, 2) == 0
        im_f(1,:,:) = 0;
        im_f(:,1,:) = 0;
    end

    im_f = reshape(im_f, [L^2 n]);

    im_f(stat_mask(:),:) = 1/2*im_f(stat_mask(:),:);
    im_f = im_f(real_mask(:),:);

    im_f = reshape(im_f, [sum(real_mask(:))*n 1]);

    vol = 1/L*anufft3(im_f, pts_rot, L*ones(1, 3));

    vol = 2*real(vol);
end
