function compute_I2E()
end

function compute_B2E()
end

function compute_In()
end

function compute_Bn()
end

function element_areas = compute_element_areas(node_data, element_data)
    % Compute the element areas from the cross product of the
    % displacement vectors
    
    % Have to do some reshaping of the view in order to do vectorized
    % indexing
    points = reshape(node_data(element_data', :), [3, size(element_data, 1), 2]);

    % Grab the two component vectors
    vec_1 = squeeze(points(2,:, :) - points(1,:, :))';
    vec_2 = squeeze(points(3,:, :) - points(1,:, :))';

    % Cross product for the area
    element_areas = 0.5 * abs(vec_1(1,:) .* vec_2(2,:) - vec_1(2,:) .* vec_2(1,:));
end


[node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("test.gri");

compute_element_areas(node_data, element_data)