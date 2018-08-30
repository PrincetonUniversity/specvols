% Helper function to print large connected components of a graph
function print_largest_components(W, min_component_size)
    [num, which] = connected_components(W);
    tab = tabulate(which);
    tab(tab(:,2) > min_component_size, 1:2)
end
