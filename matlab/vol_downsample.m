% VOL_DOWNSAMPLE Blur and downsample a volume
%
% Usage
%    vol_ds = vol_downsample(vol, N_ds, filter);
%
% Input
%    vol: Set of volumes to be downsampled in the form of an array
%       N-by-N-by-N-by-M, where M is the number of volumes.
%    N_ds: The desired resolution of the downsampled volumes.
%    filter: The type of filter to use: 'gaussian', or 'sinc' (default
%       'gaussian').
%
% Output
%    vol_ds: An array of the form N_ds-by-N_ds-by-N_ds-by-M consisting of the
%       blurred and downsampled volumes.

function vol_ds = vol_downsample(vol, N_ds, filter)
	if nargin < 3 || isempty(filter)
		filter = 'gaussian';
	end

	N = size(vol, 1);

	c = N/N_ds/(2*sqrt(2*log(2)));
	blur = exp(-[-ceil(3*c):ceil(3*c)].^2/(2*c^2));

	mesh = mesh_3d(N);
	mesh_ds = mesh_3d(N_ds);

	vol_ds = zeros([N_ds*ones(1, 3) size(vol, 4)]);

	ds_fun = @(vol)(interpn(mesh.x, mesh.y, mesh.z, vol, ...
		mesh_ds.x, mesh_ds.y, mesh_ds.z, 'linear', 0));

	if strcmp(filter, 'gaussian')
		for s = 1:size(vol_ds, 4)
			vol_s = convn(vol(:,:,:,s), blur, 'same');
			vol_ds(:,:,:,s) = ds_fun(vol_s);
		end
	else
		mask = (abs(mesh.x) < N_ds/N) & (abs(mesh.y) < N_ds/N) & ...
			(abs(mesh.z) < N_ds/N);

		vol = centered_ifft3(bsxfun(@times, centered_fft3(vol), mask));

		for s = 1:size(vol_ds, 4)
			vol_ds(:,:,:,s) = ds_fun(vol(:,:,:,s));
		end
	end
end

