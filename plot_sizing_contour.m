%creates a contour of the plot
%
%Xpts = node coordinates of a mesh
%X_BE_spline = coordinates of the splined boundary used for sizing function
function plot_sizing_contour(Xpts, X_BE_spline)

    % compute the sizing function for each node on the mesh
    for i = 1:size(Xpts,1)
        h(i) = sizing(Xpts(i,:), X_BE_spline);
    end

    %create contour plot
    x = Xpts(:, 1);
    y = Xpts(:, 2);
    [X, Y] = meshgrid(x,y);
    contour(X, Y, h);
    grid on
    title();

end