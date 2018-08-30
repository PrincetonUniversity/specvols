function test_graph_knn()
    X = [7 8;
        5 1;
        3 4;
        8 7;
        9 3];

    W = graph_knn(X,2);
    
    EXPECTED_DISTANCES_SQUARED = [
        0 0 32 2 29;
        0 0 13 0 20;
        32 13 0 0 0;
        2 0 0 0 17;
        29 20 0 17 0];
    
    assert(isequal(round(full(W).^2), EXPECTED_DISTANCES_SQUARED), "blah"); 
end