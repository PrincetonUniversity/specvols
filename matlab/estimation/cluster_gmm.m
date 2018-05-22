% CLUSTER_GMM Cluster coordinates using Gaussian mixture model
%
% Usage
%    [centers, idx, covars, taus, l_likelihood] = cluster_gmm(x, k, ...
%       cluster_opt);
%
% Input
%    x: An array of size p-by-n, containing n data points as column vectors of
%       dimension p.
%    k: The number of clusters to estimate.
%    cluster_opt: An options structure containing the fields:
%          - 'max_iter': The number of expectation-maximization iterations
%              (default 50).
%
% Output
%    centers: A p-by-k array containing the centers of the Gaussians as column
%       vectors.
%    idx: The most likely assignments of each point in one of the k clusters.
%    covars: A p-by-p-by-k array of covariance matrices for each of the
%       Gaussians.
%    taus: A vector of length k with the component probabilities of the
%       mixture.
%    l_likelihood: The log-likelihood of the parameters given the data.
%
% Note
%    The `cluster_gmm` function depends on the random number state of `rand`,
%    so to obtain reproducible results, its state must be controlled prior to
%    calling.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [centers, idx, covars, taus, l_likelihood] = cluster_gmm(x, k, ...
    cluster_opt)

    if nargin < 3 || isempty(cluster_opt)
        cluster_opt = struct();
    end

    cluster_opt = fill_struct(cluster_opt, ...
        'max_iter', 50);

    max_iter = cluster_opt.max_iter;

    p = size(x, 1);
    n = size(x, 2);

    rand('state', 1);

    taus = ones(1, k)/k;

    centers = x(:,randperm(n, k));

    covars = zeros([p*ones(1, 2) k]);
    for j = 1:k
        covars(:,:,j) = mean(var(x, [], 2))/k*eye(p);
    end

    for iter = 1:max_iter
        prob = zeros(k, n);

        for ell = 1:n
            for j = 1:k
                prob(j,ell) = gaussian_prob(x(:,ell), centers(:,j), ...
                    covars(:,:,j));
            end
        end

        l_likelihood = sum(log(sum(prob, 1)), 2);

        prob = bsxfun(@times, taus', prob);

        prob = bsxfun(@times, prob, 1./sum(prob, 1));

        prob(~isfinite(prob)) = 0;

        % what if one component has probability zero?

        taus = sum(prob, 2)'/sum(prob(:));

        for j = 1:k
            centers(:,j) = x*prob(j,:)'/sum(prob(j,:));

            c_x = bsxfun(@minus, x, centers(:,j));

            covars(:,:,j) = c_x*diag(prob(j,:))*c_x'/ ...
                sum(prob(j,:));
        end
    end

    [~, idx] = max(prob, [], 1);
end

function prob = gaussian_prob(x, center, covar)
    p = size(x, 1);

    prob = 1/((2*pi)^(p/2)*det(covar)^(1/2))* ...
        exp(-(x-center)'*inv(covar)*(x-center)/2);
end
