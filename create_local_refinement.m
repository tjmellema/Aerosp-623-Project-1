%
% Local refinement script
%

% Get the full path of the currently running script
% Setup paths to the src/ folder
pwd = fileparts(mfilename('fullpath'));
addpath(fullfile(pwd, 'src'));
% Setup pathes to the test_grids/ folder
pwd = fileparts(mfilename('fullpath'));
addpath(fullfile(pwd, 'test_grids'));

% Read in the data
[node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("coarse.gri");
%plotMesh(element_data, node_data, 'tri', true);

% Grab mappings
[I2E, B2E, In, Bn, Area] = mappings(node_data, element_data, boundary_mappings, periodic_pairs);

% Grab the coordinates of the blades
top_nodes = boundary_mappings(boundary_mappings(:,3) == 7, 1:2);
bottom_nodes = boundary_mappings(boundary_mappings(:,3) == 8, 1:2);

top_coords = unique(node_data(top_nodes,:), 'rows');
bottom_coords = unique(node_data(bottom_nodes,:), 'rows');

% Create the splines
%[spline_top, spline_bottom] = spline_boundary(top_coords, bottom_coords, 3);

[nodes, E2N] = local_refinement(node_data, element_data, top_coords, bottom_coords, I2E, B2E);
%plotMesh(E2N, nodes, 'tri', false);
