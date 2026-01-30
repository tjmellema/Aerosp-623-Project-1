%Determines the distance to all the boundary cooridnates to for projection
%
%point_cord = potential point to be added
%wall_cords = Nx2 vector of all boundary wall cords in order
%wall_distances = Nx1 vector of distances to all boundary wall cords
function wall_distances = distance2BW(point_cord, wall_cords)
    diffs = wall_cords - point_cord;
    wall_distances = sqrt(sum(diffs.^2, 2));
end