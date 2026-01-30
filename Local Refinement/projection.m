%find the point closest to the splined wall
%
%X_pt = point in the interar mesh to get rpojection for
%X_BE_spline = Boundary wall spline
function [d, xb, yb] = projection(X_pt, X_BE_spline)
    %get all distances
    distances = distance2BW(X_pt, X_BE_spline);

    %find the closest point on spline
    [d, closestIndex] = min(distances);
    xb = X_BE_spline(closestIndex,1);
    yb = X_BE_spline(closestIndex,2);
end