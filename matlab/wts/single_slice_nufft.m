%SINGLE_SLICE_NUFFT Simple wrapper function to calculate NUFFT for a single
% image at a single specified rotation
% 
% Usage
%    vol = single_slice_nufft( im, rot );
%
% Input
%    im: an L by L matrix to be "embedded" in 3 dimensions at the
%       appropriate angle
%    rot: a 3x3 matrix with determinant 1; an element of SO(3) specifying
%       the angle at which to embed the image into the volume (default
%       eye(3) )
%    fourier_flag: a binary variable that sets whether the output volume
%       will be in real space or fourier space (default no).  Currently
%       unused.
%
% Output
%   vol: vol is either an LxLxL volume, or a 2Lx2Lx2L volume, depending on
%       whether the output volume is in Fourier space or not.  It is the
%       image, embedded into 3D space at the angle specified by rot.
%       Currently always real.
%

function [ vol ] = single_slice_nufft( im, rot, fourier_flag )
    if nargin <3
        fourier_flag = 0;
    end
    if nargin < 2
        rot = eye(3)
    end
    if size(rot) ~= [3, 3]
        error('rot must be 3x3')
    end
    if abs(det(rot)) < (1-1e-4) || abs(det(rot)) > (1+1e-4)
        error('rot must have determinant 1')
    end
    
    L = size(im,1);
    
    pts_rot = rotated_grids(L,rot);
    pts_rot = reshape(pts_rot,[3 L^2]);
    im_f = 1/L^2 * centered_fft2(im);
    if mod(L, 2) == 0
        im_f(1,:,:) = 0;
        im_f(:,1,:) = 0;
    end

    im_f_flat = reshape(im_f, [L^2 1]);
    vol = 1/L*anufft3(im_f_flat, pts_rot, [L L L]);

end

