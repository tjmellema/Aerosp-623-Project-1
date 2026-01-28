import spline_curvature_ref.m
import mesh.m

%Determines the distance to all the boundary cooridnates to for projection
%
%point_cord = potential point to be added
%wall_cords = Nx2 vector of all boundary wall cords in order
%wall_distances = Nx1 vector of distances to all boundary wall cords
function wall_distances = distance2BW(point_cord, wall_cords)
    diffs = wall_cords - point_cord;
    wall_distances = sqrt(sum(diffs.^2, 2));
end

%creates universal boundary spline
%
%X_BE = original boundary points to be splined
%ref = degree of refinement
function X_BE_spline = spline_boundary(X_BE, ref)
    X_BE_spline = spline_curvature_ref.spline_curvature_ref(X_BE, length(X_BE) -1, 3, ref, 1, 1);
end

%find the point closest to the splined wall
%
%X_pt = point in the interar mesh to get rpojection for
%X_BE_spline = Boundary wall spline
function [d, xb] = projection(X_pt, X_BE_spline)
    %get all distances
    distances = distance2BW(X_pt, X_BE_spline);

    %find the closest point on spline
    [d, closestIndex] = min(distances);
    xb = X_BE_spline(closestIndex,:);

end

%find the ideal cell size at a location
%
%X_pt = point in the interar mesh to get rpojection for
%X_BE_spline = Boundary wall spline
function h = sizing(X_pt, X_BE_spline)
        %project the point on the wall
        [d, xb] = projection(X_pt, X_BE_spline);
        
        %sizing formula
        a = 1; %sixing factor
        h_max = 2; %max cell size
        % cosine spacing between leading and trailing edge
        % linear in distnace from the airfoil
        h = a.* min(h_max, (1 + cos(xb./19)) ./ 2 + d);
end


function plot_projection()
    return %TODO
end

function plot_sizing_contour()
    return %TODO
end


%Chat-GPT logic-checked and produced a cleaned version of this file
