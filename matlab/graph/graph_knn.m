% GRAPH_KNN build a symmetric knn-connectivity graph
%
% Input
%   X: n by p matrix, representing a set of n p-dimensional vectors
%   k: the k constant for the k-nearest-neighbors
%
% Output
%   W: symmetric knn-graph (n by n sparse matrix)  
%      W_{i,j} is set to 1 if if X_i is one of X_j's k nearest neighbors
%      (or vice versa), and set to 0 otherwise.
function W = graph_knn(X,k)
    [n,p] = size(X);
    [idx,D] = knnsearch(X,X,'K',k+1);
    Vi = reshape(repmat([1:n],1,k),n*k,1);
    Vj = reshape(idx(:,2:end),n*k,1);
    edge_weights = ones(length(Vi),1);
    S = sparse(Vi,Vj,edge_weights);
    W = max(S,S');
end
