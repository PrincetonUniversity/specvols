% SRC_WIENER_COORDS Calcualte coordinates using Wiener filter
%
% Usage
%    coords = src_wiener_coords(src, mean_vol, eig_vols, lambdas, noise_var);
%
% Input
%    src: A source object containing the images whose coordinates we want.
%    mean_vol: The mean volume of the source in an L-by-L-by-L array.
%    eig_vols: The eigenvolumes of the source in an L-by-L-by-L-by-K array.
%    lambdas: The eigenvalues in a K-by-K diagonal matrix (default `eye(K)`).
%    noise_var: The variance of the noise in the images (default 0).
%
% Output
%    coords: A K-by-`src.n` array of coordinates corresponding to the Wiener
%       filter coordinates of each image in `src`. These are obtained by
%       formula
%
%          alpha_s = eig_vols^T H_s ( y_s - P_s mean_vol ) ,
%
%       where P_s is the forward image mapping and y_s is the sth image,
%
%          H_s = Sigma * P_s^T ( P_s Sigma P_s^T + noise_var I )^(-1) ,
%
%       and Sigma is the covariance matrix eig_vols * lambdas * eig_vols^T.
%       Note that when noise_var is zero, this reduces to the projecting
%       y_s onto the span of P_s eig_vols.

% TODO: Enable batch processing.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function coords = src_wiener_coords(src, mean_vol, eig_vols, lambdas, ...
    noise_var)

    if nargin < 5 || isempty(noise_var)
        noise_var = 0;
    end

    if nargin < 4 || isempty(lambdas)
        lambdas = eye(size(eig_vols, 4));
    end

    [Qs, Rs] = qr_vols_forward(src, 1, src.n, eig_vols);

    ims = src_image(src, 1, src.n);

    ims = ims - vol_forward(src, mean_vol, 1, src.n);

    covar_noise = noise_var*eye(size(eig_vols, 4));

    coords = zeros([size(eig_vols, 4), src.n], src.precision);

    for s = 1:src.n
        im_coords = im_to_vec(Qs(:,:,:,s))'*im_to_vec(ims(:,:,s));

        covar_im = (Rs(:,:,s)*lambdas*Rs(:,:,s)' + covar_noise);

        im_coords = lambdas*Rs(:,:,s)'*(covar_im\im_coords);

        coords(:,s) = im_coords;
    end
end

function [Qs, Rs] = qr_vols_forward(src, s, n, vols)
    ims = zeros([src.L*ones(1, 2) n size(vols, 4)], class(vols));

    for ell = 1:size(vols, 4)
        ims(:,:,:,ell) = vol_forward(src, vols(:,:,:,ell), s, n);
    end

    ims = permute(ims, [1 2 4 3]);

    Qs = zeros([src.L*ones(1, 2) size(vols, 4)], class(vols));
    Rs = zeros([size(vols, 4)*ones(1, 2) n], class(vols));

    for k = 1:n
        [Qk, Rk] = qr(im_to_vec(ims(:,:,:,k)), 0);

        Qs(:,:,:,k) = vec_to_im(Qk);
        Rs(:,:,k) = Rk;
    end
end
