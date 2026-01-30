%
% Create the coarse mesh
%

% Some parameters from the project spec
inlet_height = 18;
inlet_length  = 17;

% Simplex element scaling. Adjust this to get ~1000–1500 triangles.
% This is multipled by the chord length for the maximum element size.
h_max_scale = 1/15;

% Tolerance factor for identical points 
tol = 1e-10;

% Blade file names
upper_blade_file = 'contrib/bladeupper.txt';
lower_blade_file = 'contrib/bladelower.txt';

% Load the data and run some assert checks on it
assert(isfile(upper_blade_file), "Can't find %s", upper_blade_file);
assert(isfile(lower_blade_file), "Can't find %s", lower_blade_file);
upper_blade_data = load(upper_blade_file);
lower_blade_data = load(lower_blade_file);
assert(size(upper_blade_data, 2) == 2 && size(lower_blade_data, 2) == 2, 'Blade files must be nx2 arrays');

% Shift the lower surface up by the inlet height
lower_blade_data(:, 2) = lower_blade_data(:, 2) + inlet_height;

% For creating the other lines in the mesh we need to know the min and
% max points of the blade data. Thus, we have to sort the data based 
% on x position. 
% NOTE: TBH we don't technically have to sort the data, we can just take
% the argmin and argmax but I'll be lazy.
upper_blade_data = sortrows(upper_blade_data, 1);
lower_blade_data = sortrows(lower_blade_data, 1);

min_lower_blade = lower_blade_data(1, :);
min_upper_blade = upper_blade_data(1, :);
max_lower_blade = lower_blade_data(end, :);
max_upper_blade = upper_blade_data(end, :);

% We can check that the points match for the x positions
assert(min_lower_blade(1) == min_upper_blade(1), "min x points don't match for blade data");
assert(max_lower_blade(1) == max_upper_blade(1), "max x points don't match for blade data");

% Compute the chord length for the element scale later on
chord_length = max_lower_blade(1) - min_lower_blade(1);

% Now compute some points for the inlet and outlet construction
inlet_top = min_lower_blade - [inlet_length, 0];
inlet_bottom = min_upper_blade - [inlet_length, 0];
outlet_top = max_lower_blade + [inlet_length, 0];
outlet_bottom = max_upper_blade + [inlet_length, 0];

% Now we assemble the surface loop. Importantly, order matters for this
% one. As such, we'll traverse counter-clockwise and flip the lower blade
lower_blade_data = flipud(lower_blade_data);

surface_loop = [ upper_blade_data; ...
                 outlet_bottom; ...
                 outlet_top; ...
                 lower_blade_data; ...
                 inlet_top; ...
                 inlet_bottom; ...
               ];

% Plot the surface loop 
% Note this isn't actually closed due to not wanting duplicate
% points
if do_plots
    figure('Name','Surface loop'); hold on; axis equal; grid on;
    plot(surface_loop(:,1), surface_loop(:,2), 'b.-');
    xlabel('x (mm)'); ylabel('y (mm)');
end

% Create the pde model and geometry objects
pde_model = createpde();
surface_loop_x = surface_loop(:,1);
surface_loop_y = surface_loop(:,2);
geometry_description = [2; numel(surface_loop_x); surface_loop_x(:); surface_loop_y(:)];

geometry_matrix = decsg(geometry_description);
geometry = geometryFromEdges(pde_model, geometry_matrix);

if do_plots
    figure('Name','Edge labels');
    pdegplot(pde_model, 'EdgeLabels','on'); axis equal; grid on;
end

% Generate the coarse mesh from the objects above
% Compute the max element size
h_max = chord_length * h_max_scale;
msh = generateMesh(pde_model, 'Hmax', h_max, 'GeometricOrder', 'linear');

if do_plots
    figure('Name','Coarse mesh');
    pdemesh(pde_model); axis equal; grid on;
    title(sprintf('Coarse mesh: Hmax=%.4g, Triangles=%d', h_max, size(msh.Elements, 2)));
end

fprintf('Hmax = %.6g\n', h_max);
fprintf('Triangles = %d\n', size(msh.Elements,2));
