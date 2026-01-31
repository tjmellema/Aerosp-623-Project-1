function h = sizing(point, top_spline, bottom_spline)
    %
    % Sizing function for a measure of how big to mesh cells
    %

    % First determine the distance from the spline of interest
    [top_distance, top_point] = projection(point, top_spline);
    [bottom_distance, bottom_point] = projection(point, bottom_spline);

    % Determine which spline the point ought to belong to
    if top_distance < bottom_distance
        distance = top_distance; 
        x_coord = top_point(:, 1); 
    else
        distance = bottom_distance; 
        x_coord = bottom_point(:, 1); 
    end

        
    % Sizing formula
    a = .545; %sizing factor
    h_max = 4; %max cell size
    h_min = .01; %min cell size

    % cosine spacing between leading and trailing edge
    % linear in distance from the airfoil
    h = a.* max(min(h_max, (1 + cos(x_coord./19)) ./ 2 + .25*distance), h_min);
end