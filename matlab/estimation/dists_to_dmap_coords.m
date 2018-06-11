% DISTS_TO_DMAP_COORDS Calculate diffusion map coordinates from distances
%
% Usage
%    dmap_coords = dists_to_dmap_coords(dists, epsilon, t, num_coords);
%
% Input
%    dists: An n-by-n positive matrix containing distances.
%    espilon: The kernel width (default `median(dists(:).^2)`).
%    num_coords: The number of diffusion map coordinates to calculate
%       (default 3).
%    t: The diffusion time (default 0).
%
% Output
%    dmap_coords: The calculated diffusion map coordinates in an array of size
%       num_coords-by-n.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function dmap_coords = dists_to_dmap_coords(dists, epsilon, num_coords, t)
    if nargin < 4 || isempty(t)
        t = 0;
    end

    if nargin < 3 || isempty(num_coords)
        num_coords = 3;
    end

    if nargin < 2 || isempty(epsilon)
        epsilon = median(dists(:).^2);
    end

    W = exp(-dists.^2/epsilon);

    D = diag(sum(W, 2));

    S = D^(-1/2)*W*D^(-1/2);

    S = make_symmat(S);

    [V, lambda] = eigs(S, num_coords+1, 'la');

    phi = D^(-1/2)*V;

    dmap_coords = phi(:,2:end)*lambda(2:end,2:end)^t;
    dmap_coords = dmap_coords.';
end
