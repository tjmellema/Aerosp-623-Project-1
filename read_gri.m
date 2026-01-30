function [node_data, element_data, boundary_mappings, periodic_pairs] = read_gri(filename)
    % File stream
    fstream = fopen(filename, 'r');

    % Read in global data from the first three points
    A = fscanf(fstream,'%d', 3);
    n_nodes = A(1);
    %n_element_ids = A(2);
    dim = A(3);
    assert(dim == 2, "Only 2D is supported")

    % Allocate an array to store the point data for the nodes
    node_data = zeros(n_nodes, dim);

    % Read in node data
    for i = 1:n_nodes
        point = fscanf(fstream, '%lf', dim);
        node_data(i, :) = point;
    end

    % Read in the number of unique boundary groups (ids).
    A = fscanf(fstream,'%d', 1);
    n_boundary_ids = A(1);

    % Read the boundary info
    boundary_mappings = zeros(0, 3);
    for i = 1:n_boundary_ids
        % Grab the overall data for the boundary group (id)
        A = textscan(fstream,'%d %d %s', 1);
        n_boundary_face = A{1};
        %n_linear_nodes = A{2};
        %boundary_type = A{3};

        % Grab the rest of the nodes that are part of the boundary
        for j = 1:n_boundary_face
            boundary_face_nodes = fscanf(fstream, '%d', 2);
            boundary_mappings(end+1, :) = [boundary_face_nodes.', i];
        end
    end
    
    % Read in the element group info
    element_data = zeros(0, 3);
    % TODO: Fix this element group stuff so it isn't hard-coded
    % in. Shouldn't matter much for this project though
    n_element_groups = 1;
    for i = 1:n_element_groups
        % Read in the data for the element group
        A = textscan(fstream,'%d %d %s', 1);
        n_elements = A{1};
        %element_order = A{2};
        basis_function = A{3};

        assert(strcmp(basis_function, 'TriLagrange'), "Invalid basis function")
        % Read in the element nodes for each element
        for j = 1:n_elements
            element_nodes = fscanf(fstream, '%d', 3);
            element_data(end+1, :) = element_nodes;
        end
    end

    % Read in the periodicity data
    periodic_pairs = zeros(0, 2);
    A = textscan(fstream,'%d %s', 1);
    n_periodic_groups = A{1};
    assert(strcmp(A{2}, 'PeriodicGroup'), "Invalid .gri syntax for periodicity")
    for i = 1:n_periodic_groups
        A = textscan(fstream,'%d %s', 1);
        n_periodic_node_pairs = A{1};
        periodicity_type = A{2};
        assert(strcmp(periodicity_type,'Translational'), "Rotational periodicity not implemented")

        for j = 1:n_periodic_node_pairs
            periodic_pair = fscanf(fstream, '%d', 2);
            periodic_pairs(end+1, :) = periodic_pair;
        end
    end
end