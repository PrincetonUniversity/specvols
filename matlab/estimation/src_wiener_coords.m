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
%    coords_opt: An options structure containing the fields:
%          - 'batch_size': The size of the batches in which to compute the
%             coordinates (default 512).
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

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function coords = src_wiener_coords(src, mean_vol, eig_vols, lambdas, ...
    noise_var, coords_opt)

    if nargin < 6 || isempty(coords_opt)
        coords_opt = struct();
    end

    if nargin < 5 || isempty(noise_var)
        noise_var = 0;
    end

    if nargin < 4 || isempty(lambdas)
        lambdas = eye(size(eig_vols, 4));
    end

    coords_opt = fill_struct(coords_opt, ...
        'batch_size', 512);

    coords = zeros([size(eig_vols, 4), src.n], src.precision);

    covar_noise = noise_var*eye(size(eig_vols, 4));

    for batch = 1:ceil(src.n/coords_opt.batch_size)
        batch_s = (batch-1)*coords_opt.batch_size+1;
        batch_n = min(batch*coords_opt.batch_size, src.n)-batch_s+1;

        [Qs, Rs] = qr_vols_forward(src, batch_s, batch_n, eig_vols);

        ims = src_image(src, batch_s, batch_n);

        ims = ims - vol_forward(src, mean_vol, batch_s, batch_n);

        for s = 1:batch_n
            im_coords = im_to_vec(Qs(:,:,:,s))'*im_to_vec(ims(:,:,s));

            covar_im = (Rs(:,:,s)*lambdas*Rs(:,:,s)' + covar_noise);

            im_coords = lambdas*Rs(:,:,s)'*(covar_im\im_coords);

            coords(:,batch_s+s-1) = im_coords;
        end
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
