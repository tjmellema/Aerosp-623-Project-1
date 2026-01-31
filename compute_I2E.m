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