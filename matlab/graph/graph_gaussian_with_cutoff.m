% GRAPH_GAUSSIAN_WITH_CUTOFF build a sparse weighted graph based on an
%   almost gaussian kernel with compact support
%
% Input
%   X: n by p matrix, representing p-dimensional vectors
%   sigma: standard deviation of the gaussian kernel
%   num_sigma_cutoff: disconnect pairs of vertices that are at a distance
%   greater than sigma*num_sigma_cutoff
%
% Output:
%   W: sparse graph weights satisfying (n by n sparse matrix)
%
%                 { e^(-||X_i-X_j||^2 / 2 sigma^2) if ||X_i-X_j|| < sigma*num_sigma_cutoff
%       W_{i,j} = { 
%                 { 0                                 otherwise
%
% Description
%   The purpose of this construction is to be close to the gaussian kernel construction
%   implemented in graph_gaussian_kernel while producing a sparse graph.
%
%   This implementation is a combination of graph_epsilon and graph_gaussian_kernel.

function W = graph_gaussian_with_cutoff(X, sigma, num_sigma_cutoff)
    [n, p] = size(X);
    [idx,dist] = rangesearch(X, X, num_sigma_cutoff*sigma);
    row_index_as_cell = num2cell(1:length(idx));
    Vi = cell2mat(cellfun(@(row, i) i*ones(1,length(row)), idx', row_index_as_cell, 'UniformOutput', false));
    Vj = cell2mat(idx');
    distances_ij = cell2mat(dist');
    weights = exp(-distances_ij.^2 / (2*sigma^2));
    W = sparse(Vi, Vj, weights);
end
 
