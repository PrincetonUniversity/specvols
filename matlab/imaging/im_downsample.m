% IM_DOWNSAMPLE Blur and downsample image
%
% Usage
%    im_ds = im_downsample(im, L_ds);
%
% Input
%    im: Set of images to be downsampled in the form of an array L-by-L-by-K,
%       where K is the number of images.
%    L_ds: The desired resolution of the downsampled images. Must be smaller
%       than L.
%
% Output
%    im_ds: An array of the form L_ds-by-L_ds-by-K consisting of the blurred
%       and downsampled images.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function im_ds = im_downsample(im, L_ds)
    N = size(im, 1);

    grid = grid_2d(N);
    grid_ds = grid_2d(L_ds);

    im_ds = zeros([L_ds*ones(1, 2) size(im, 3)], class(im));

    % NOTE: Cannot use interp2 because it's not compatible with ndgrid.
    ds_fun = @(im)(interpn(grid.x, grid.y, im, grid_ds.x, grid_ds.y, ...
        'linear', 0));

    mask = (abs(grid.x) < L_ds/N) & (abs(grid.y) < L_ds/N);

    im = real(centered_ifft2(bsxfun(@times, centered_fft2(im), mask)));

    for s = 1:size(im_ds, 3)
        im_ds(:,:,s) = ds_fun(im(:,:,s));
    end
end
