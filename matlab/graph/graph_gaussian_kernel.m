% GRAPH_GAUSSIAN_KERNEL build a weighted graph based on a gaussian kernel
%
% Input
%   X: n by p matrix, representing p-dimensional vectors
%   sigma: standard deviation of the gaussian kernel
%
% Output
%   W: graph weights satisfying W_{i,j} = e^(-||X_i-X_j||^2 / 2 sigma^2)
%      (n by n full matrix)

function W = graph_gaussian_kernel(X, sigma)
    pairwise_distances = squareform(pdist(X));
    W = exp(-pairwise_distances.^2 / (2*sigma^2));
    for i=1:length(W)
        W(i,i) = 0;
    end
end
