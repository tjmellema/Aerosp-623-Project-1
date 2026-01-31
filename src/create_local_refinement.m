%
% Local refinement script
%

% Read in the data
[node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("coarse.gri");

% Grab the coordinates of the blades
top_nodes = boundary_mappings(boundary_mappings(:,3) == 7, 1:2);
bottom_nodes = boundary_mappings(boundary_mappings(:,3) == 8, 1:2);

top_coords = unique(node_data(top_nodes,:), 'rows');
bottom_coords = unique(node_data(bottom_nodes,:), 'rows');

% Create the splines
[spline_top, spline_bottom] = spline_boundary(top_coords, bottom_coords, 3);

% Build the 
[IE, BE] = BuildEdges(element_data);

[nodes, E2N] = localRefinement(node_data,element_data,BE,IE, spline_top, spline_bottom, periodic_pairs);
