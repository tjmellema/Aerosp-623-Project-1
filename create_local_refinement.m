%local refinement call

[node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("passage_coarse.gri");

top_nodes = unique(boundary_mappings(boundary_mappings(:,3)==7,1:2));
bot_nodes = unique(boundary_mappings(boundary_mappings(:,3)==8,1:2));

top_coords = node_data(top_nodes,:);
bot_coords = node_data(bot_nodes,:);

[X_BE_spline_top, X_BE_spline_bottom] = spline_boundary(top_coords, bot_coords, 3);
[IE, BE] = BuildEdges(element_data);
[nodes, E2N] = localRefinement(node_data,element_data,BE,IE, X_BE_spline_top, X_BE_spline_bottom);
