% VOL_DOWNSAMPLE Blur and downsample volume
%
% Usage
%    vol_ds = vol_downsample(vol, L_ds);
%
% Input
%    vol: Set of volumes to be downsampled in the form of an array
%       L-by-L-by-L-by-K, where K is the number of volumes.
%    L_ds: The desired resolution of the downsampled volumes. Must be smaller
%       than L.
%
% Output
%    vol_ds: An array of the form L_ds-by-L_ds-by-L_ds-by-K consisting of the
%       blurred and downsampled volumes.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vol_ds = vol_downsample(vol, L_ds)
    N = size(vol, 1);

    grid = grid_3d(N);
    grid_ds = grid_3d(L_ds);

    vol_ds = zeros([L_ds*ones(1, 3) size(vol, 4)], class(vol));

    % NOTE: Cannot use interp2 because it's not compatible with ndgrid.
    ds_fun = @(vol)(interpn(grid.x, grid.y, grid.z, vol, ...
        grid_ds.x, grid_ds.y, grid_ds.z, 'linear', 0));

    mask = (abs(grid.x) < L_ds/N) & (abs(grid.y) < L_ds/N) & (abs(grid.z) < L_ds/N);

    vol = real(centered_ifft3(bsxfun(@times, centered_fft3(vol), mask)));

    for s = 1:size(vol_ds, 4)
        vol_ds(:,:,:,s) = ds_fun(vol(:,:,:,s));
    end
end
