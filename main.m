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
write_gri("passage_coarse.gri", nodes, elements, ...
          inlet_height, inlet_length, ...
          inlet_top, inlet_bottom, ...
          outlet_top, outlet_bottom ...
         );

% Verify meshes
mesh_verification('test_grids/test.gri');
mesh_verification('test_grids/test_2.gri');
mesh_verification('test_grids/test_3.gri');
mesh_verification('passage_coarse.gri');

create_local_refinement();

write_gri("passage_coarse_local_refinement.gri", nodes, E2N, ...
          inlet_height, inlet_length, ...
          inlet_top, inlet_bottom, ...
          outlet_top, outlet_bottom ...
         );
%{

mesh_verification('test_grids/test1.gri');

mesh_verification('test_grids/test2.gri');
mesh_verification('test_grids/testModified.gri');
mesh_verification('test_grids/testModified1.gri');
mesh_verification('test_grids/testModified2s.gri');
%}