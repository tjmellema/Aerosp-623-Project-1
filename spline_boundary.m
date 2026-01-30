%creates universal boundary spline
%
%X_BE = original boundary points to be splined
%ref = degree of refinement
function X_BE_spline = spline_boundary(X_BE, ref)
    X_BE_spline = spline_curvature_ref(X_BE, length(X_BE) -1, 3, ref, 1, 1);
end