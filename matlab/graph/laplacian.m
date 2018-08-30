% LAPLACIAN compute the graph laplacian operator given a graph's weights matrix.
%
% Input
%   W: graph weights, given as a symmetric sparse matrix.
%      You can use any of the graph_* functions for building such matrices.
%   laplacian_type: one of the following: 'combinatorial', 'randomwalk', 'normalized', 'geometric'
%
% Output
%   L: sparse matrix containing the graph Laplacian
%
% Description
%   Let D the degree matrix, a diagonal matrix satisfying D_ii = \sum_j W_ij
%   This function computes one of the following Laplacians:
%   'combinatorial'
%       L = D - W
%   'randomwalk'
%       Lw = D^{-1} L = I - D^{-1} W
%   'normalized'
%       L_norm = D^{-1/2} L D^{-1/2} = I - D^{-1/2} W D^{-1/2} 
%   'geometric'
%       L = laplacian(D^{-1} L D^{-1}, 'randomwalk')
%
% Notes
%   1. If the data points are sampled on a manifold, the geometric Laplacian
%      advocated by Coifman & Lafon (2006) is the only one that converges to the
%      Laplace-Beltrami operator, regardless of the density.
%   2. Before computing the laplacian, we zero out the diagonal of W.
%      This has been shown by El Karui & Hu (2016) to make the laplacian more robust to noise.
%
% References
%   [1] Coifman & Lafon (2006) "Diffusion maps"
%   [2] El Karoui and Wu (2016) "Graph connection laplacian methods can be made robust to noise"

function L = laplacian(W, laplacian_type)
    [nrows, ncols] = size(W);
    assert(nrows == ncols);
    n = nrows;
    
    % Zero out diagonal
    W(1:1+n:end) = 0;

    degrees = sum(W,2);

    switch laplacian_type
        case 'combinatorial'
            D = spdiags(degrees, [0], n, n);
            L = D-W;

        case 'randomwalk'
            I = speye(n);
            assert(all(degrees > 0), 'Vertex degrees must all be positive');
            L = I - bsxfun(@rdivide, W, degrees'); % normalize W's rows

        case 'normalized'
            I = speye(n);
            assert(all(degrees > 0), 'Vertex degrees must all be positive');
            sqrtdegrees = sqrt(degrees);
            L = I - bsxfun(@rdivide, bsxfun(@rdivide, W, sqrtdegrees), sqrtdegrees');
            
        case 'geometric'
            assert(all(degrees > 0), 'Vertex degrees must all be positive');
            %W_new = (W./degrees)./degrees';
            W_new = bsxfun(@rdivide, bsxfun(@rdivide, W, degrees), degrees');
            L = laplacian(W_new, 'randomwalk');

        otherwise
            error(message('laplacian:invalid laplacian type'));
    end
end

