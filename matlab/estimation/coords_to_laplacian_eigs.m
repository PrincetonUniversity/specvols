% COORDS_TO_LAPLACIAN_EIGS Calculate eigenvectors of graph laplacian from
% distance matrix
%
% Usage
%    laplacian_eigs = coords_to_laplacian_eigs(dists, epsilon, num_coords, t);
%
% Input
%    coords: An k-by-n positive matrix containing coordinates of points on 
%    manifold.
%    espilon: The kernel width (default `median(dists(:).^2)`).
%    num_coords: The number of diffusion map coordinates to calculate
%       (default 3).
%    t: The diffusion time (default 0).
%
% Output
%    laplacian_eigs: The calculated diffusion map coordinates in an array of size
%       num_coords-by-n.

% Author
%    Joakim Anden <janden@flatironinstitute.org>
%    Amit Halevi <ahalevi@princeton.edu>

function laplacian_eigs = coords_to_laplacian_eigs(coords, epsilon, num_coords)

    if nargin < 3 || isempty(num_coords)
        num_coords = 3;
    end

   
    n2 = sum(abs(coords).^2, 1);

    d2 = bsxfun(@plus, n2, n2') - 2*coords'*coords;

    clear n2
    
    d2_thresholded = max(0,d2);
    
    clear d2;

    if nargin < 2 || isempty(epsilon)
        epsilon = median(d2_thresholded(:)) / 4;
    end
  
    dists = sqrt(d2_thresholded);
    
    clear d2_thresholded
  
    sq_dists = dists.^2;
    
    clear dists
    
    norm_dists = sq_dists/(-epsilon);
    
    clear sq_dists
    
    W = exp(norm_dists);

    clear norm_dists
    
    D = diag(sum(W, 2));

    d = diag(D);

    clear D
        
    S_not_sym = bsxfun(@times,1./d,W);

    clear W
    
    S = make_symmat(S_not_sym);
    
    clear S_not_sym
    
    [V, lambda] = eigs(double(S), num_coords, 'la');

    phi = bsxfun(@times,1./sqrt(d),V);

    laplacian_eigs = phi';
end
