import spline_curvature_ref.m

%point_cord = potential point to be added
%wall_cords = Nx2 vector of all boundary wall cords in order
%wall_distances = Nx1 vector of distances to all boundary wall cords

%Chat-GPT logic-checked and produced a cleaned version of this function

%Determines the distance to all the boundary cooridnates to for projection
function wall_distances = distance2BW(point_cord, wall_cords)
    diffs = wall_cords - point_cord;
    wall_distances = sqrt(sum(diffs.^2, 2));
end

%TODO - get splined boundary BEFORE computing the projection so we can
%reuse the spline in each iteration.
function [d, xb] = projection(point_cord, wall_cords)
    %get all distances
    distances = distance2BW(point_cord, wall_cords);

    %find the closest point
    [minDistance, closestIndex] = min(distances);
    
    %find the closest  between left and right point
    rightBool = 1;
    if distances(closestIndex-1) < distances(closestIndex+1)
        rightBool = 0;
    end

    %spline between pt1 and pt2
    Xnew = spline_curvature_ref.spline_curvature_ref(wallcords((closestIndex-2+rightBool-8):(cloesestIndex+1+rightBool+8)), ...
                                                     20, 3, 3, 1, 1);

    %find the closest distance along spline
    spline_distances = distance2BW(point_cord, Xnew);

    %return d and xb from spline
    [d, d_i] = min(spline_distances);
    xb = Xnew(d_i, 1);

end

function h = sizing(sizing_point, wall_points)
        %project the point on the wall
        [d, xb] = projection(ptIN, wall_points);
        
        %sizing formula
        h = d + xb; %TODO
end