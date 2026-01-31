function [distance, point] = projection(point, spline)
    %
    % Find the point closest to a spline
    %
    distances = distance_to_spline(point, spline);
    [distance, closest_index] = min(distances);
    point = spline(closest_index, :);
end