function test_graph_epsilon()
    X = rand(100,3);
    epsilon = 1;
    distances = squareform(pdist(X));
    W = spgraph_epsilon(X, epsilon);
    assert(isequal(W, distances < epsilon));
end
