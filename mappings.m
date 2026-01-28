function face_pairs = compute_face_pairs(element_data)
    %    
    % Compute all the faces pairs of the given node and element data
    %

    % Extract face pairs from the element data
    face_pairs = [element_data(:, [1, 2]); element_data(:, [2, 3]); element_data(:, [3, 1])];
    face_pairs = unique(sort(face_pairs, 2), 'rows');
end

function mask = compute_periodic_face_pair_mask(face_pairs, periodic_pairs)
    %
    % Compute the bool mask for the periodic face pairs
    %

    % For the periodic pairs, we have the node mappings. However, we
    % need to identify the faces from that. To do so, we'll do the 
    % forward substitution of the mappings to identity the identical pairs.
    
    % First grab the number of nodes
    n_nodes = max(face_pairs(:));
    % Construct the map (identity for now)
    forward_map = (1:n_nodes).';
    % Update the map with the periodic pairs 1 -> 2
    forward_map(periodic_pairs(:, 1)) = periodic_pairs(:, 2);
    % Substitute
    faces_forward_sub = forward_map(face_pairs);

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

function interior_face_to_edge = compute_I2E(node_data, element_data, boundary_mappings, periodic_pairs)
    % When constructing the interior face data, we need to fill in all
    % but the boundary mappings, except for the periodic face pairs.
    
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
    sorted_interior_faces = sort(interior_faces, 2)
    sorted_element_data = sort(element_data, 2)

    % Find the elements that match the interior edge

end

function compute_B2E()
end

function compute_In()
end

function compute_Bn()
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

plotgri("test.gri")
[node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("test.gri");
compute_I2E(node_data, element_data, boundary_mappings, periodic_pairs)