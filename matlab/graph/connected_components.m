% CONNECTED_COMPONENTS Split a graph into its connected components
%
% Input
%   G: nxn sparse matrix representing the graph weights
%
% Output
%   n_conected_components:
%       Number of connected components
%   which:
%       integer vector of length n such that which(i) is the
%       connected component number of the i-th vertex.
%
% Description
%   Uses the built-in Dulmage-Mendelsohn decomposition which (as a side effect)
%   decomposes the graph into blocks.

function [n_connected_component, which] = connected_components(G)
    [p,q,block_indices] = dmperm(G+speye(size(G)));
    n_connected_component = numel(block_indices)-1;
    block_boundary_indicator = full(sparse(1,block_indices,1));
    block_counter = cumsum(block_boundary_indicator);
    which(p) = block_counter(1:end-1);
end
