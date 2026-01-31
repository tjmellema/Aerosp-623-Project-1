function distance = distance_to_spline(point_cord, spline)
    %
    % Find distance to a spline
    %
    displacement = spline - point_cord;
    distance = sqrt(sum(displacement.^2, 2));
end