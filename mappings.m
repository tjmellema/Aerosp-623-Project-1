function [I2E, B2E, In, Bn, Area] = mappings(node_data, element_data, boundary_mappings, periodic_pairs)
    % wrapper function
    
    % 1. Compute Connectivity Matrices
    I2E = compute_I2E(element_data, boundary_mappings, periodic_pairs);
    I2E = unique(I2E, "rows");
    B2E = compute_B2E(element_data, boundary_mappings, periodic_pairs);
    
    % 2. Compute Normals
    In = compute_In(node_data, element_data, I2E);
    Bn = compute_Bn(node_data, element_data, B2E);
    
    % 3. Compute Areas
    Area = compute_element_areas(node_data, element_data);
end

function face_pairs = compute_face_pairs(element_data)
    %    
    % Compute all the faces pairs of the given node and element data
    %

    % Extract face pairs from the element data
    face_pairs = [element_data(:, [1, 2]); element_data(:, [2, 3]); element_data(:, [3, 1])];
    face_pairs = unique(sort(face_pairs, 2), 'rows');
end

function data_forward_sub = compute_forward_substitution(data, periodic_pairs)
    %
    % Compute the forward substitution of the data given some periodic
    % mappings.
    %

     % First grab the number of nodes
    n_nodes = max(data(:));
    % Construct the map (identity for now)
    forward_map = (1:n_nodes).';
    % Update the map with the periodic pairs 1 -> 2
    forward_map(periodic_pairs(:, 1)) = periodic_pairs(:, 2);
    % Substitute
    data_forward_sub = forward_map(data);
end

function data_backward_sub = compute_backward_substitution(data, periodic_pairs)
    %
    % Compute the backward substitution of the data given some periodic
    % mappings.
    %

     % First grab the number of nodes
    n_nodes = max(data(:));
    % Construct the map (identity for now)
    forward_map = (1:n_nodes).';
    % Update the map with the periodic pairs 2 -> 1
    forward_map(periodic_pairs(:, 2)) = periodic_pairs(:, 1);
    % Substitute
    data_backward_sub = forward_map(data);
end

function mask = compute_periodic_face_pair_mask(face_pairs, periodic_pairs)
    %
    % Compute the bool mask for the periodic face pairs
    %

    % For the periodic pairs, we have the node mappings. However, we
    % need to identify the faces from that. To do so, we'll do the 
    % forward substitution of the mappings to identity the identical pairs.
    faces_forward_sub = compute_forward_substitution(face_pairs, periodic_pairs);

    % Now find the mask of identical faces.
    % NOTE: Order matters here! Given the previous substitution, we get
    % some degenerate faces that can be misidentified if we don't care 
    % about order.
    [~, ~, unique_in_array_indices] = unique(faces_forward_sub, 'rows', 'stable');
    pair_occurence_count = accumarray(unique_in_array_indices, 1);
    mask = pair_occurence_count(unique_in_array_indices) > 1;
end

function boundary_face_pairs = compute_boundary_face_pairs(boundary_mappings, periodic_pairs)
    %    
    % Compute all the boundary fair pairs (excluding periodic faces)
    %

    % All the boundary face pairs are given here
    face_pairs = boundary_mappings(:, 1:2);

    % Periodic face pair mask
    periodic_face_mask = compute_periodic_face_pair_mask(face_pairs, periodic_pairs);

    % Slice the non-periodic ones
    boundary_face_pairs = face_pairs(~periodic_face_mask, :);
end

function mask = compute_node_mask(element_nodes, face_nodes)
    %
    % Compute the node mask of the given face nodes on the element nodes
    %
    face_node_1 = face_nodes(:, 1);
    face_node_2 = face_nodes(:, 2);
    mask = (element_nodes == face_node_1) | (element_nodes == face_node_2);
end

function interior_face_to_edge = compute_I2E(element_data, boundary_mappings, periodic_pairs)
    % When constructing the interior face data, we need to fill in all
    % but the boundary mappings, except for the periodic face pairs.
    % NOTE: We technically have a degenerate state for the matched interior
    % edge for the periodic boundary condition. This can be solved with a
    % unqiue call I believe. Hopefully it's not an issue right now.
    
    % Grab the fair pairs
    faces = compute_face_pairs(element_data);

    % Grab the boundary face pairs
    % NOTE: This ignores periodic boundary pairs
    boundary_faces =  compute_boundary_face_pairs(boundary_mappings, periodic_pairs);
    
    % Boundary face mask
    is_boundary_face = ismember(faces, boundary_faces, 'rows');

    % Get only the interior face node pairs
    interior_faces = faces(~is_boundary_face, :);

    % Preallocate the array now
    interior_face_to_edge = zeros(size(interior_faces, 1), 4);

    % Sort the interior face and element data arrays
    sorted_interior_faces = sort(interior_faces, 2);
    sorted_element_data = sort(element_data, 2);

    % Substitute
    element_data_forward_sub = compute_forward_substitution(sorted_element_data, periodic_pairs);
    % Re-sort 
    sorted_element_data_forward_sub = sort(element_data_forward_sub, 2);

    % Substitute
    element_data_backward_sub = compute_backward_substitution(sorted_element_data, periodic_pairs);
    % Re-sort 
    sorted_element_data_backward_sub = sort(element_data_backward_sub, 2);

    % Find the elements that match the interior edge
    % NOTE: Hopefully this loop isn't to slow...
    for i = 1:size(sorted_interior_faces, 1)
        node_pair = sorted_interior_faces(i, :);
        node_pair_forward_sub = compute_forward_substitution(node_pair, periodic_pairs)';
        node_pair_backward_sub = compute_backward_substitution(node_pair, periodic_pairs)';

        % Check which elements correspond to this face pair. We should 
        % have two elements, so if that criterion isn't met check for 
        % periodic mappings
        edge_1 = ismember(sorted_element_data(:, [1,2]), node_pair, 'rows');
        edge_2 = ismember(sorted_element_data(:, [1,3]), node_pair, 'rows');
        edge_3 = ismember(sorted_element_data(:, [2,3]), node_pair, 'rows');
        total_edges = edge_1 + edge_2 + edge_3;
        total_found_elements = sum(total_edges);

        if total_found_elements == 2
            element_indices = find(total_edges);
            interior_face_to_edge(i, [1, 3]) = element_indices;
        elseif total_found_elements == 1
            % When we have periodic pairs, we need to do a forward mapping
            % substitution on the element data, resort it, and recompute
            % masks for finding the indices. We must also do the same
            % for the backward substitution

            % Mask
            edge_1_forward = ismember(sorted_element_data_forward_sub(:, [1,2]), node_pair, 'rows');
            edge_2_forward = ismember(sorted_element_data_forward_sub(:, [1,3]), node_pair, 'rows');
            edge_3_forward = ismember(sorted_element_data_forward_sub(:, [2,3]), node_pair, 'rows');
            edge_1_backward = ismember(sorted_element_data_backward_sub(:, [1,2]), node_pair, 'rows');
            edge_2_backward = ismember(sorted_element_data_backward_sub(:, [1,3]), node_pair, 'rows');
            edge_3_backward = ismember(sorted_element_data_backward_sub(:, [2,3]), node_pair, 'rows');
            total_edges = edge_1_forward + edge_2_forward + edge_3_forward + edge_1_backward + edge_2_backward + edge_3_backward;

            element_indices = find(total_edges);
            interior_face_to_edge(i, [1, 3]) = element_indices;
        else
            assert(false, "I don't know how you got here. 2 periodic faces!")
        end


        % Now that we have the element ids we can find the corresponding 
        % local face id of the interior edge. Importantly, since the local
        % face is opposite the local node, we can find the id of the node 
        % left out.
        left_element_nodes = element_data(interior_face_to_edge(i, 1), :);
        right_element_nodes = element_data(interior_face_to_edge(i, 3), :);

        left_mask = compute_node_mask(left_element_nodes, node_pair);
        right_mask = compute_node_mask(right_element_nodes, node_pair);

        left_local_indices = find(~left_mask);
        right_local_indices = find(~right_mask);

        % Deal with periodic boundary conditions
        % TODO: I'm not entirely sure that I should only consider the 
        % backward substitution for the left element and the opposite
        % for the right element. It seems to work, but I can't guarantee 
        % it.
        if size(left_local_indices, 2) == 1
            interior_face_to_edge(i, 2) = left_local_indices;
        else
            left_mask = compute_node_mask(left_element_nodes, node_pair_backward_sub);
            left_local_indices = find(~left_mask);
            interior_face_to_edge(i, 2) = left_local_indices;
        end
        if size(right_local_indices, 2) == 1
            interior_face_to_edge(i, 4) = right_local_indices;
        else
            right_mask = compute_node_mask(right_element_nodes, node_pair_forward_sub);
            right_local_indices = find(~right_mask);
            interior_face_to_edge(i, 4) = right_local_indices;
        end
    end
end

function boundary_face_to_edge = compute_B2E(element_data, boundary_mappings, periodic_pairs)
    % Grab the boundary face pairs
    % NOTE: This ignores periodic boundary pairs
    boundary_faces =  compute_boundary_face_pairs(boundary_mappings, periodic_pairs);

    % Allocate the array
    boundary_face_to_edge = zeros(size(boundary_faces, 1), 3);
    
    % Since we don't have to worry about periodic faces, the rest becomes
    % rather simple compared to compute_I2E. First, find the corresponding
    % element id

    % Sort the face and element data arrays
    sorted_boundary_faces = sort(boundary_faces, 2);
    sorted_element_data = sort(element_data, 2);
    sorted_boundary_mappings(:, 1:2) = sort(boundary_mappings(:, 1:2), 2);

    for i = 1:size(boundary_faces, 1)
        boundary_face = sorted_boundary_faces(i, :);
        % Compute element face masks
        edge_1 = ismember(sorted_element_data(:, [1,2]), boundary_face, 'rows');
        edge_2 = ismember(sorted_element_data(:, [1,3]), boundary_face, 'rows');
        edge_3 = ismember(sorted_element_data(:, [2,3]), boundary_face, 'rows');
        total_edges = edge_1 + edge_2 + edge_3;
    
        % Element indices
        element_indices = find(total_edges);
    
        % Fill in the indices
        boundary_face_to_edge(i, 1) = element_indices;

        % Find the local node mask
        local_node_mask = compute_node_mask(element_data(element_indices, :), boundary_face);
        boundary_face_to_edge(i, 2) = find(~local_node_mask);

        % Fill in the boundary id
        mask = sum(sorted_boundary_mappings(:, 1:2) == boundary_face, 2) == 2;
        boundary_face_to_edge(i, 3) = boundary_mappings(mask, 3);
    end
end

function interior_normals = compute_In(node_data, element_data, I2E)
    %
    % Compute the normal vector going from the left to right element
    % for interior faces
    %

    % Preallocate the array
    interior_normals = zeros(size(I2E, 1), 2);

    % Grab some slices of things
    left_elements = I2E(:, 1);
    right_elements = I2E(:, 3);

    left_nodes = element_data(left_elements,  :);
    right_nodes = element_data(right_elements,  :);

    % I want to grab the positions but matrices...
    left_nodes_flat = left_nodes(:);
    right_nodes_flat = right_nodes(:);

    left_positions_flat = node_data(left_nodes_flat, :);
    right_positions_flat = node_data(right_nodes_flat, :);

    left_node_positions = reshape(left_positions_flat, [size(left_nodes, 1), 3, 2]);
    right_node_positions = reshape(right_positions_flat, [size(right_nodes, 1), 3, 2]);

    % Compute the centroid of the left and right elements
    left_node_centroid = squeeze(mean(left_node_positions, 2));
    right_node_centroid = squeeze(mean(right_node_positions, 2));

    % Displacement
    displacement = right_node_centroid - left_node_centroid;

    % Norm
    interior_normals(:, :) = displacement ./ vecnorm(displacement, 2, 2);
end

function boundary_normals = compute_Bn(node_data, element_data, B2E)
    %
    % Compute the normal vector going from the left to right element
    % for boundary faces
    %

    % Preallocate the array
    boundary_normals = zeros(size(B2E, 1), 2);

    % Grab some slices of things
    elements = B2E(:, 1);
    faces = B2E(:, 2);
    nodes = element_data(elements,  :);
    face_nodes = nodes(sub2ind(size(nodes), (1:size(nodes,1))', faces));
    boundary_face_nodes = nodes(~(nodes == face_nodes));


    % I want to grab the positions but matrices...
    boundary_positions_flat = node_data(boundary_face_nodes, :);
    boundary_positions = reshape(boundary_positions_flat, [2, size(nodes, 1), 2]);

    % Displacement
    displacement = squeeze(boundary_positions(2, :, :) - boundary_positions(1, :, :));

    % Norm
    boundary_normals(:, [2,1]) = displacement ./ vecnorm(displacement, 2, 2);

end

function element_areas = compute_element_areas(node_data, element_data)
    %
    % Compute the element areas from the cross product of the
    % displacement vectors
    %

    % Have to do some reshaping of the view in order to do vectorized
    % indexing
    points = reshape(node_data(element_data', :), [3, size(element_data, 1), 2]);

    % Grab the two component vectors
    vec_1 = squeeze(points(2,:, :) - points(1,:, :))';
    vec_2 = squeeze(points(3,:, :) - points(1,:, :))';

    % Cross product for the area
    element_areas = 0.5 * abs(vec_1(1,:) .* vec_2(2,:) - vec_1(2,:) .* vec_2(1,:));
end

% plotgri("test.gri")
% [node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("test.gri");
% I2E = compute_I2E(element_data, boundary_mappings, periodic_pairs);
% B2E = compute_B2E(element_data, boundary_mappings, periodic_pairs);
% interior_normals = compute_In(node_data, element_data, I2E);
% boundary_normals = compute_Bn(node_data, element_data, B2E);