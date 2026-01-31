function write_grid(filename, nodes, elements, ...
                               inlet_height, inlet_length, ...
                               inlet_top, inlet_bottom, ...
                               outlet_top, outlet_bottom ...
                              )
    n_nodes = size(nodes, 1);
    n_elements = size(elements, 1);
    dim = size(nodes, 2);

    assert(dim == 2, "Unsupported dimension")

    % Create a triangulation object given the nodes and elements
    tria = triangulation(elements, nodes);

    % Grab the boundary faces
    boundary_faces = freeBoundary(tria);

    % Grab the boundary face points
    % The first dimension is the face pair, the second is all the points,
    % and the last dimension is the spatial dimension.
    boundary_faces_flat = boundary_faces(:);
    boundary_face_points = reshape(nodes(boundary_faces_flat, :), [], dim, dim);

    % Grab the boundary faces from each point set mask. We'll start with
    % the periodic one since that's the most important.
    tol = 1e-10;

    x = boundary_face_points(:,:,1);
    y = boundary_face_points(:,:,2);

    on_inlet_top = abs(y - inlet_top(2)) < tol;
    on_inlet_bottom = abs(y - inlet_bottom(2)) < tol;
    on_outlet_top = abs(y - outlet_top(2)) < tol;
    on_outlet_bottom = abs(y - outlet_bottom(2)) < tol;

    on_inlet_side = abs(x - inlet_top(1)) < tol;
    on_outlet_side = abs(x - outlet_top(1)) < tol;

    left_of_blade = x <= inlet_top(1) + inlet_length + tol;
    right_of_blade = x >= outlet_top(1) - inlet_length - tol;

    inlet_top_mask = sum(on_inlet_top & left_of_blade, 2) == 2;
    inlet_bottom_mask = sum(on_inlet_bottom & left_of_blade, 2) == 2;
    outlet_top_mask = sum(on_outlet_top & right_of_blade, 2) == 2;
    outlet_bottom_mask = sum(on_outlet_bottom & right_of_blade, 2) == 2;

    inlet_side_mask = sum(on_inlet_side & left_of_blade, 2) == 2;
    outlet_side_mask = sum(on_outlet_side & right_of_blade, 2) == 2;

    combined_mask = inlet_top_mask | inlet_bottom_mask | ...
                    outlet_top_mask | outlet_bottom_mask | ...
                    inlet_side_mask | outlet_side_mask;

    blade_boundary_faces = boundary_faces(~combined_mask, :);
    blade_boundary_points = boundary_face_points(~combined_mask, :, :);

    x = blade_boundary_points(:,:,1);
    y = blade_boundary_points(:,:,2);


    % Janky masking because I don't know an elegant way to do it
    cond1 = (x < 3) & (y > 10);
    cond2 = (x >= 2) & (x < 6) & (y > 7);
    cond3 = (x >= 4) & (y > 3);

    blade_top_mask = ...
      (sum(cond1, 2) == 2) | ...
      (sum(cond2, 2) == 2) | ...
      (sum(cond3, 2) == 2);
    blade_bottom_mask = ~blade_top_mask;

    % Finding the periodic node matches. Because we are in 2D this is as 
    % simple as matching node pairs with the same x coordinate.
    inlet_top_nodes = unique(boundary_faces(inlet_top_mask, :));
    inlet_bottom_nodes = unique(boundary_faces(inlet_bottom_mask, :));
    assert(size(inlet_top_nodes, 1) == size(inlet_bottom_nodes, 1), "Invalid periodicity node matching");

    inlet_top_xy    = nodes(inlet_top_nodes, :);
    inlet_bottom_xy = nodes(inlet_bottom_nodes, :);
    
    outlet_top_nodes = unique(boundary_faces(outlet_top_mask, :));
    outlet_bottom_nodes = unique(boundary_faces(outlet_bottom_mask, :));
    assert(size(outlet_top_nodes, 1) == size(outlet_bottom_nodes, 1), "Invalid periodicity node matching");

    outlet_top_xy    = nodes(outlet_top_nodes, :);
    outlet_bottom_xy = nodes(outlet_bottom_nodes, :);

    [~, index_inlet_top] = sort(inlet_top_xy(:,1));
    [~, index_inlet_bottom] = sort(inlet_bottom_xy(:,1));
    
    [~, index_outlet_top] = sort(outlet_top_xy(:,1));
    [~, index_outlet_bottom] = sort(outlet_bottom_xy(:,1));
    
    inlet_top_nodes    = inlet_top_nodes(index_inlet_top);
    inlet_bottom_nodes = inlet_bottom_nodes(index_inlet_bottom);
    
    outlet_top_nodes    = outlet_top_nodes(index_outlet_top);
    outlet_bottom_nodes = outlet_bottom_nodes(index_outlet_bottom);

    % Check that the x positions of the node match 1 to 1. If not, 
    % shift them.
    tol = 1e-12;
    x_top    = nodes(inlet_top_nodes, 1);
    x_bottom = nodes(inlet_bottom_nodes, 1);
    diff_x = abs(x_top - x_bottom);
    if ~all(diff_x < tol)
        warning('Periodic nodes mismatch in x-coordinates! Shifting them');

        % Compute the average x for each pair
        x_avg = (x_top + x_bottom) / 2;
        % Replace values
        nodes(inlet_top_nodes, 1)    = x_avg;
        nodes(inlet_bottom_nodes, 1) = x_avg;
    end
    x_top    = nodes(outlet_top_nodes, 1);
    x_bottom = nodes(outlet_bottom_nodes, 1);
    diff_x = abs(x_top - x_bottom);
    if ~all(diff_x < tol)
        warning('Periodic nodes mismatch in x-coordinates! Shifting them');

        % Compute the average x for each pair
        x_avg = (x_top + x_bottom) / 2;
        % Replace values
        nodes(outlet_top_nodes, 1)    = x_avg;
        nodes(outlet_bottom_nodes, 1) = x_avg;
    end

    % Open an output stream
    fstream = fopen(filename, 'w');

    % Write the header info
    fprintf(fstream, '%d %d %d\n', n_nodes, n_elements, dim);

    % Write the node data
    for i = 1:n_nodes
        fprintf(fstream, '%g ', nodes(i, :));
        fprintf(fstream, '\n');
    end

    % Write the number of unique boundary groups (ids).
    % With our given geometry we have 8 groups
    n_boundary_ids = 8;
    n_linear_nodes = dim;
    fprintf(fstream, '%d\n', n_boundary_ids);
    for i = 1:n_boundary_ids
        switch i
            case 1
                n_boundary_face = sum(inlet_top_mask);
                local_boundary_faces = boundary_faces(inlet_top_mask, :);
                boundary_type = "InletTop";
            case 2
                % NOTE: We have to flip the direction here
                n_boundary_face = sum(inlet_bottom_mask);
                local_boundary_faces = flip(boundary_faces(inlet_bottom_mask, :), 2);
                boundary_type = "InletBottom";
            case 3
                n_boundary_face = sum(outlet_top_mask);
                local_boundary_faces = boundary_faces(outlet_top_mask, :);
                boundary_type = "OutletTop";
            case 4
                % NOTE: We have to flip the direction here
                n_boundary_face = sum(outlet_bottom_mask);
                local_boundary_faces = flip(boundary_faces(outlet_bottom_mask, :), 2);
                boundary_type = "OutletBottom";
            case 5
                n_boundary_face = sum(inlet_side_mask);
                local_boundary_faces = boundary_faces(inlet_side_mask, :);
                boundary_type = "InletSide";
            case 6
                n_boundary_face = sum(outlet_side_mask);
                local_boundary_faces = boundary_faces(outlet_side_mask, :);
                boundary_type = "OutletSide";
            case 7
                n_boundary_face = sum(blade_top_mask);  
                local_boundary_faces = blade_boundary_faces(blade_top_mask, :);
                boundary_type = "BladeTop";
            case 8
                n_boundary_face = sum(blade_bottom_mask);  
                local_boundary_faces = blade_boundary_faces(blade_bottom_mask, :);
                boundary_type = "BladeBottom";
            otherwise
                assert(false, "How did you get here?");
        end

        fprintf(fstream, '%d %d %s\n', n_boundary_face, n_linear_nodes, boundary_type);

        for j = 1:n_boundary_face
            fprintf(fstream, '%d %d\n', local_boundary_faces(j, :));
        end
    end

    % Write the element group info
    % Since MATLAB does all the hard work of order the elements with 
    % counter-clockwise connectivity, we don't have to do anything special
    % here. We also don't have anything special with the elements because
    % everything is trilgrange and first order.
    n_element_groups = 1;
    element_order = 1;
    basis_function = 'TriLagrange';

    for i = 1:n_element_groups
        fprintf(fstream, '%d %d %s\n', n_elements, element_order, basis_function);
        for j = 1:n_elements
            fprintf(fstream, '%d %d %d\n', elements(j,1), elements(j,2), elements(j,3));
        end
    end

    % Write the periodic mapping info
    n_periodic_groups = 2;
    
    periodicity_type = 'Translational';
    fprintf(fstream, '%d PeriodicGroup\n', n_periodic_groups);
    for i = 1:n_periodic_groups
        switch i
            case 1
                n_periodic_node_pairs = size(inlet_top_nodes, 1);
                mappings = [inlet_top_nodes inlet_bottom_nodes];
            case 2
                n_periodic_node_pairs = size(outlet_top_nodes, 1);
                mappings = [outlet_top_nodes outlet_bottom_nodes];
            otherwise
                assert(false, "How did you get here?");
        end
        fprintf(fstream, '%d %s\n', n_periodic_node_pairs, periodicity_type);
            for j = 1:n_periodic_node_pairs
                fprintf(fstream, '%d %d\n', mappings(j, 1), mappings(j, 2));
            end
    end
    
    % Clean up
    fclose(fstream);
end