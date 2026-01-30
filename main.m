%
% Main function that runs each subpart of the code.
%

clear;
close all;
clc;

% Setup pathes to the src/ folder
pwd = fileparts(mfilename('fullpath'));
addpath(fullfile(pwd, 'src'));

% Generate the coarse grid
do_plots = false;
mesh();

% Write coarse grid to file
write_gri("passage_coarse.gri", msh, ...
          inlet_height, inlet_length, ...
          inlet_top, inlet_bottom, ...
          outlet_top, outlet_bottom ...
         );

plotgri("passage_coarse.gri")

% Verify meshes
mesh_verification('passage_coarse.gri');
%{
mesh_verification('test_grids/test.gri');
mesh_verification('test_grids/test_2.gri');

mesh_verification('test_grids/test1.gri');

mesh_verification('test_grids/test2.gri');
mesh_verification('test_grids/testModified.gri');
mesh_verification('test_grids/testModified1.gri');
mesh_verification('test_grids/testModified2s.gri');
%}