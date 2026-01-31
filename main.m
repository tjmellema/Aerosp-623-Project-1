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
mesh_verification('test_grids/test.gri');
mesh_verification('test_grids/test_2.gri');
mesh_verification('test_grids/test_3.gri');
mesh_verification('coarse.gri');


create_local_refinement();

write_gri("coarse_local_refinement.gri", nodes, E2N, ...
          inlet_height, inlet_length, ...
          inlet_top, inlet_bottom, ...
          outlet_top, outlet_bottom ...
         );