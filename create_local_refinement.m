%
% Local refinement script
%

% Read in the data
[node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("coarse.gri");

top_nodes = boundary_mappings(boundary_mappings(:,3) == 7, 1:2);
bottom_nodes = boundary_mappings(boundary_mappings(:,3) == 8, 1:2);

top_coords = unique(node_data(top_nodes,:), 'rows');
bot_coords = unique(node_data(bottom_nodes,:), 'rows');

[X_BE_spline_top, X_BE_spline_bottom] = spline_boundary(top_coords, bot_coords, 3);
[IE, BE] = BuildEdges(element_data);
[nodes, E2N] = localRefinement(node_data,element_data,BE,IE, X_BE_spline_top, X_BE_spline_bottom, periodic_pairs);
