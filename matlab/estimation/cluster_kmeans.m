% CLUSTER_KMEANS Perform k-means clustering using Lloyd's algorithm
%
% Usage
%    [centers, idx] = cluster_kmeans(x, k, cluster_opt);
%
% Input
%    x: An array of size p-by-n, containing n data points as column vectors of
%       dimension p.
%    k: The desired number of clusters.
%    cluster_opt: An options structure containing the fields:
%          - 'max_iter': The maximum number of iterations to run for each
%             clustering attempt (default 128).
%          - 'replicates': The number of clustering attempts to perform, after
%             which the clustering with the lowest mean square distance is
%             selected (default 8).
%
% Output
%    centers: A p-by-k array containing the centers of the clusters.
%    idx: An 1-by-n array of indices corresponding to the cluster assigned to
%       each data point.
%
% Note
%    The `cluster_kmeans` function depends on the random number state of
%    `rand`, so to obtain reproducible results, its state must be controlled
%    prior to calling.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [centers, idx] = cluster_kmeans(x, k, cluster_opt)
    if nargin < 3
        cluster_opt = struct();
    end

    cluster_opt = fill_struct(cluster_opt, ...
        'max_iter', 128, ...
        'replicates', 8);

    x2 = sum(abs(x).^2, 1);

    d2_best = Inf;
    centers_best = zeros(size(x, 1), k);
    idx_best = zeros(size(x,1), 2);

    for rep = 1:cluster_opt.replicates
        centers = x(:,randperm(size(x, 2), k));

        for iter = 1:cluster_opt.max_iter
            centers2 = sum(abs(centers).^2, 1);

            d2 = bsxfun(@plus, x2, centers2') - 2*centers'*x;

            [d2, idx] = min(d2, [], 1);

            d2 = sum(d2);

            for ell = 1:k
                centers(:,ell) = mean(x(:,idx==ell), 2);
            end
        end

        if d2 < d2_best
            centers_best = centers;
            idx_best = idx;

            d2_best = d2;
        end
    end

    centers = centers_best;
    idx = idx_best;
end
