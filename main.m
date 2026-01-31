%
% Main function that runs each subpart of the code.
%

clear;
close all;
clc;

% Setup paths to the src/ folder
pwd = fileparts(mfilename('fullpath'));
addpath(fullfile(pwd, 'src'));
% Setup pathes to the Local Refinement/ folder
pwd = fileparts(mfilename('fullpath'));
addpath(fullfile(pwd, 'Local Refinement'));


% Generate the coarse grid
do_plots = false;
mesh();

% Grab some data from the FEMesh object
nodes = msh.Nodes.';
elements = msh.Elements.';

% Write coarse grid to file
write_gri("coarse.gri", nodes, elements, ...
          inlet_height, inlet_length, ...
          inlet_top, inlet_bottom, ...
          outlet_top, outlet_bottom ...
         );

% Verify meshes
tic;
mesh_verification('test_grids/test.gri');
time1 = toc;
fprintf('Time for test.gri: %.4f seconds\n', time1);

tic;
mesh_verification('test_grids/test_2.gri');
time2 = toc;
fprintf('Time for test_2.gri: %.4f seconds\n', time2);

tic;
mesh_verification('test_grids/test_3.gri');
time3 = toc;
fprintf('Time for test_3.gri: %.4f seconds\n', time3);

tic;
mesh_verification('coarse.gri');
time4 = toc;
fprintf('Time for coarse.gri: %.4f seconds\n', time4);

tic;
mesh_verification('passage_coarse1.gri');
time5 = toc;
fprintf('Time for passage_coarse1.gri: %.4f seconds\n', time5);

tic;
mesh_verification('passage_coarse2.gri');
time6 = toc;
fprintf('Time for passage_coarse2.gri: %.4f seconds\n', time6);

tic;
mesh_verification('passage_coarse3.gri');
time7 = toc;
fprintf('Time for passage_coarse3.gri: %.4f seconds\n', time7);


create_local_refinement();
plotMesh(E2N, nodes, 'tri')

write_gri("coarse_local_refinement.gri", nodes, E2N, ...
          inlet_height, inlet_length, ...
          inlet_top, inlet_bottom, ...
          outlet_top, outlet_bottom ...
         );
