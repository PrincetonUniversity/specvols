% GAUSSIAN_BLOB_VOLS Generate Gaussian blob volumes
%
% Usage
%    vols = gaussian_blob_vols(L, C, K, alpha);
%
% Input:
%    L: The size of the volumes (default 8).
%    C: The number of volumes to generate (default 2).
%    K: The number of blobs (default 16).
%    alpha: A scale factor of the blob widths (default 1).
%    precision: The precision of the volumes, 'single' or 'double' (default
%       'single').
%
% Output
%    vols: A volume array of size L-by-L-by-L-by-C containing the C Gaussian
%       blob volumes.
%
% Note
%    The 'gaussian_blob_vols' function depends on the random number state of
%    'randn' so to obtain reproducible results, its state must be controlled
%    prior to calling.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vols = gaussian_blob_vols(L, C, K, alpha, precision)
    if nargin < 1 || isempty(L)
        L = 8;
    end

    if nargin < 2 || isempty(C)
        C = 2;
    end

    if nargin < 3 || isempty(K)
        K = 16;
    end

    if nargin < 4 || isempty(alpha)
        alpha = 1;
    end

    if nargin < 5 || isempty(precision)
        precision = 'single';
    end

    vols = zeros([L*ones(1, 3) C], precision);

    for k = 1:C
        [Q, D, mu] = gaussian_blobs(K, alpha, precision);

        vols(:,:,:,k) = eval_gaussian_blobs(L, Q, D, mu);
    end
end

function [Q, D, mu] = gaussian_blobs(K, alpha, precision)
    Q = zeros([3*ones(1, 2) K], precision);
    D = zeros([3*ones(1, 2) K], precision);
    mu = zeros([3 K], precision);

    for k = 1:K
        V = randn(3, 3, precision)/sqrt(3);
        [Q(:,:,k), ~] = qr(V);
        D(:,:,k) = alpha^2/16*diag(sum(abs(V).^2, 1));

        mu(:,k) = 0.5*randn(3, 1, precision)/sqrt(3);
    end
end

function vol = eval_gaussian_blobs(L, Q, D, mu)
    grid3d = grid_3d(L);

    coords = [grid3d.x(:) grid3d.y(:) grid3d.z(:)]';

    K = size(Q, 3);

    vol = zeros([1 size(coords, 2)], class(Q));

    for k = 1:K
        coords_k = bsxfun(@minus, coords, mu(:,k));
        coords_k = Q(:,:,k)*1/sqrt(D(:,:,k))*Q(:,:,k)'*coords_k;

        vol = vol+exp(-sum(abs(coords_k).^2, 1)/2);
    end

    vol = reshape(vol, size(grid3d.x));
end
