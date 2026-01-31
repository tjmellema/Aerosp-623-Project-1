% Setup paths to the src/ folder
pwd = fileparts(mfilename('fullpath'));
addpath(fullfile(pwd, 'src'));
addpath(fullfile(pwd, 'test_grids'));

[node_data, E2N, boundary_mappings, periodic_pairs] = read_gri('passage_coarse1.gri');
[I2E, B2E, In, Bn, Area] = mappings(node_data, E2N, boundary_mappings, periodic_pairs);

                %(Xs, num_iter, E2N, B2E, num_nodes)
Xsnew = Smooth(node_data, 5, E2N, B2E, size(node_data, 1));