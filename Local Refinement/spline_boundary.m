%creates universal boundary spline
%
%X_BE = original boundary points to be splined
%ref = degree of refinement
function [X_BE_spline_top, X_BE_spline_bottom] = spline_boundary(X_BE_top, X_BE_bottom, ref)
    X_BE_spline_top = spline_curvature_ref(X_BE_top, length(X_BE_top) -1, 3, ref, 1, 1);
    X_BE_spline_bottom = spline_curvature_ref(X_BE_bottom, length(X_BE_bottom) -1, 3, ref, 1, 1);
end