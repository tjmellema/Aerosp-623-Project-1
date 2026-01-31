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
    % Grab the face pairs
    faces = compute_face_pairs(element_data);
    boundary_faces = sort(compute_boundary_face_pairs(boundary_mappings, periodic_pairs), 2);
    
    % Boundary face mask
    is_boundary_face = ismember(faces, boundary_faces, 'rows');
    assert(size(boundary_faces, 1) == sum(is_boundary_face), "Mismatch in boundary face size")
    
    interior_faces = faces(~is_boundary_face, :);
    interior_face_to_edge = zeros(size(interior_faces, 1), 4);
    
    sorted_interior_faces = sort(interior_faces, 2);
    sorted_element_data = sort(element_data, 2);
    
    % Pre-compute substitutions
    elem_fwd = sort(compute_forward_substitution(sorted_element_data, periodic_pairs), 2);
    elem_back = sort(compute_backward_substitution(sorted_element_data, periodic_pairs), 2);
    % --- OPTIMIZATION: Linear Indexing (Pairing) ---
    % Convert N x 2 rows into a single unique integer: hash = min*K + max
    % This replaces 'ismember(..., 'rows')' with a simple vector comparison.
    max_node = max([element_data(:); sorted_interior_faces(:)]);
    K = max_node + 1; 
    
    % Hash the element edges
    e1 = sorted_element_data(:,1)*K + sorted_element_data(:,2);
    e2 = sorted_element_data(:,1)*K + sorted_element_data(:,3);
    e3 = sorted_element_data(:,2)*K + sorted_element_data(:,3);
    
    % Hash the periodic element edges
    ef1 = elem_fwd(:,1)*K + elem_fwd(:,2); ef2 = elem_fwd(:,1)*K + elem_fwd(:,3); ef3 = elem_fwd(:,2)*K + elem_fwd(:,3);
    eb1 = elem_back(:,1)*K + elem_back(:,2); eb2 = elem_back(:,1)*K + elem_back(:,3); eb3 = elem_back(:,2)*K + elem_back(:,3);
    for i = 1:size(sorted_interior_faces, 1)
        node_pair = sorted_interior_faces(i, :);
        % Hash the target node pair
        target_hash = node_pair(1)*K + node_pair(2);
        
        % Check standard edges
        total_edges = (e1 == target_hash) + (e2 == target_hash) + (e3 == target_hash);
        found_count = sum(total_edges);
        
        if found_count == 2
            element_indices = find(total_edges);
            interior_face_to_edge(i, [1, 3]) = element_indices;
        else
            % Check periodic edges (Forward and Backward)
            total_edges_p = (ef1 == target_hash) + (ef2 == target_hash) + (ef3 == target_hash) + ...
                            (eb1 == target_hash) + (eb2 == target_hash) + (eb3 == target_hash);
            
            % Combine standard and periodic matches
            all_matches = total_edges | total_edges_p;
            element_indices = find(all_matches);
            
            % Sanity check: Should always have exactly 2 elements for an interior/periodic face
            if length(element_indices) ~= 2
                 error("Face %d connects to %d elements. Check periodic connectivity.", i, length(element_indices));
            end
            interior_face_to_edge(i, [1, 3]) = element_indices;
        end
        % Node masking logic (exactly as original)
        node_pair_fwd = compute_forward_substitution(node_pair, periodic_pairs)';
        node_pair_back = compute_backward_substitution(node_pair, periodic_pairs)';
        
        left_nodes = element_data(interior_face_to_edge(i, 1), :);
        right_nodes = element_data(interior_face_to_edge(i, 3), :);
        
        % Local index logic... (rest of original code follows)
        l_mask = compute_node_mask(left_nodes, node_pair);
        r_mask = compute_node_mask(right_nodes, node_pair);
        l_idx = find(~l_mask); r_idx = find(~r_mask);
        if length(l_idx) == 1, interior_face_to_edge(i, 2) = l_idx;
        else
            l_mask = compute_node_mask(left_nodes, node_pair_back);
            l_idx = find(~l_mask);
            if length(l_idx) == 1, interior_face_to_edge(i, 2) = l_idx;
            else
                l_mask = compute_node_mask(left_nodes, node_pair_fwd);
                interior_face_to_edge(i, 2) = find(~l_mask);
            end
        end
        if length(r_idx) == 1, interior_face_to_edge(i, 4) = r_idx;
        else
            r_mask = compute_node_mask(right_nodes, node_pair_fwd);
            r_idx = find(~r_mask);
            if length(r_idx) == 1, interior_face_to_edge(i, 4) = r_idx;
            else
                r_mask = compute_node_mask(right_nodes, node_pair_back);
                interior_face_to_edge(i, 4) = find(~r_mask);
            end
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
    % Grab some slices of things

    % Note to get the right direction of the interior normals, we have to 
    % rely on the fact that the node positions are ordered counter-clock-
    % wise. Using this, we need to take the displacement at the 
    % preceding node minus the next one.
    left_elements = I2E(:, 1);
    left_faces = I2E(:, 2);
    left_nodes = element_data(left_elements,  :);
    left_face_nodes = left_nodes(sub2ind(size(left_nodes), (1:size(left_nodes,1))', left_faces));
    mask = left_nodes == left_face_nodes;
    interior_face_nodes = zeros(size(left_nodes,1), 2);
    for i = 1:size(left_nodes,1)
        zero_indices = find(~mask(i, :));
        % Edge case for where 0 is in the middle
        if zero_indices(2) - zero_indices(1) > 1
            zero_indices([2, 1]) = zero_indices([1, 2]);
        end
        interior_face_nodes(i,:) = left_nodes(i, zero_indices);
    end 
    
    % I want to grab the positions but matrices...
    interior_positions_flat = node_data(interior_face_nodes, :);
    interior_positions = reshape(interior_positions_flat, [size(left_nodes, 1), 2, 2]);

    % Displacement
    displacement = squeeze(interior_positions(:, 2, :) - interior_positions(:, 1, :));

    % Orthogonalize by rotation 90 degree clockwise
    displacement = [displacement(:,2), -displacement(:,1)];

    % Norm
    interior_normals = displacement ./ vecnorm(displacement, 2, 2);
end

function boundary_normals = compute_Bn(node_data, element_data, B2E)
    %
    % Compute the normal vector going from the left to right element
    % for boundary faces.
    %

    % Grab some slices of things

    % Note to get the right direction of the boundary normals, we have to 
    % rely on the fact that the node positions are ordered counter-clock-
    % wise. Using this, we need to take the displacement at the 
    % preceding node minus the next one.
    elements = B2E(:, 1);
    faces = B2E(:, 2);
    nodes = element_data(elements,  :);
    face_nodes = nodes(sub2ind(size(nodes), (1:size(nodes,1))', faces));
    mask = nodes == face_nodes;
    boundary_face_nodes = zeros(size(nodes,1), 2);
    for i = 1:size(nodes,1)
        zero_indices = find(~mask(i, :));
        % Edge case for where 0 is in the middle
        if zero_indices(2) - zero_indices(1) > 1
            zero_indices([2, 1]) = zero_indices([1, 2]);
        end
        boundary_face_nodes(i,:) = nodes(i, zero_indices);
    end
    
    % I want to grab the positions but matrices...
    boundary_positions_flat = node_data(boundary_face_nodes, :);
    boundary_positions = reshape(boundary_positions_flat, [size(nodes, 1), 2, 2]);

    % Displacement
    displacement = squeeze(boundary_positions(:, 2, :) - boundary_positions(:, 1, :));
    
    % Orthogonalize
    displacement = [displacement(:,2), -displacement(:,1)];

    % Norm
    boundary_normals = displacement ./ vecnorm(displacement, 2, 2);
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
    signed_area = 0.5 * ( ...
        vec_1(1,:) .* vec_2(2,:) - ...
        vec_1(2,:) .* vec_2(1,:) );

    assert(all(signed_area > 0), "Connectivities not counter-clockwise");

    element_areas = abs(signed_area);
end
