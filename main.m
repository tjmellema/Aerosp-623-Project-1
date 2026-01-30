%
% Main function that runs each subpart of the code.
%

% Setup pathes to the src/ folder
pwd = fileparts(mfilename('fullpath'));
addpath(fullfile(pwd, 'src'));

% Generate the coarse grid
mesh();

% Write coarse grid to file
WriteGriFixed();

% Verify meshes
mesh_verification('test_grids/test.gri');
mesh_verification('test_grids/test_2.gri');
mesh_verification('passage_coarse.gri');
mesh_verification('test_grids/test1.gri');
%{
mesh_verification('test_grids/test2.gri');
mesh_verification('test_grids/testModified.gri');
mesh_verification('test_grids/testModified1.gri');
mesh_verification('test_grids/testModified2s.gri');
%}