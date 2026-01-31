function h = sizing(point, top_spline, bottom_spline)
    %
    % Sizing function for a measure of how big to mesh cells
    %

    % First determine the distance from the spline of interest
    [top_distance, top_point] = projection(point, top_spline);
    [bottom_distance, bottom_point] = projection(point, bottom_spline);

    % Determine which spline the point ought to belong to
    if top_distance < bottom_distance
        d = top_distance; 
        xb = top_point(:, 1); 
    else
        d = bottom_distance; 
        xb = bottom_point(:, 1); 
    end

        
    % Sizing formula
    a = 0.545; %sizing factor
    h_max = 1000000; %max cell size
    h_min = 0; %min cell size

    % cosine spacing between leading and trailing edge
    % linear in distance from the airfoil
    % edge refinement
    A = 0.7; xLE = -9.5; xTE = 9.5; sigma = 1;
    edge_factor = 1 ...
        - A * exp(-((xb - xLE).^2) / sigma^2) ...
        - A * exp(-((xb - xTE).^2) / sigma^2);

    edge_factor = max(edge_factor, 0.3);

    h = a.* max(min(h_max, edge_factor + 0.1 .* d), h_min);
end