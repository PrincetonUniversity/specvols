% GRAPH_EPSILON build an epsilon-connectivity graph
%
% Input
%   X: n by p matrix, representing p-dimensional vectors
%   epsilon: the distance threshold. Any 2 row vectors from X of distance less than epsilon
%       will be connected.
%
% Output
%   W: epsilon-graph adjacency matrix (n by n sparse matrix)

function W = graph_epsilon(X, epsilon)
    [n, p] = size(X);
    idx = rangesearch(X, X, epsilon);

    % The following fancy cell functions are a faster way to do the following:
    % W = sparse(n,n);
    % for i = 1:length(idx)
    %     idx_row = idx{i};
    %     for j = 1:length(idx_row)
    %         W(i,idx_row(j)) = 1;
    %     end
    % end
    row_index_as_cell = num2cell(1:length(idx));
    Vi = cell2mat(cellfun(@(row, i) i*ones(1,length(row)), idx', row_index_as_cell, 'UniformOutput', false));
    Vj = cell2mat(idx');
    W = sparse(Vi, Vj, ones(1,length(Vi)));
end
 
