%plot the sizing function
[node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("test.gri");
%TODO - get X_BE!!!
X_BE_spline = spline_boundary(X_BE, 2);
plot_sizing_contour(node_data, X_BE_spline);


%requested function in 1.4 to show the projection function
%
%Xpts = set of coordinates to project
%X_BE = unsplined boundary coordinates
function plot_projection(Xpts, X_BE, ref)

    %boundary spline
    X_BE_spline = spline_boundary(X_BE, ref);

    %projection
    [d, xb, yb] = projection(Xpts, X_BE_spline);

    %plot the reference spline
    plot(X_BE_spline, Linewidth=1);
    hold on

    %plot the points
    scatter(Xpts,'-b' , Linewidth=3);

    %plot the projection to the spline
    X_plot = [xb'; X_BE_spline(1)'];
    Y_plot = [yb'; X_BE_spline(2)'];
    plot(X_plot, Y_plot, '-b', Linewidth=2);

    title('Sample Projection Function')
    grid on


end

%creates a contour of the plot
%
%Xpts = node coordinates of a mesh
%X_BE_spline = coordinates of the splined boundary used for sizing function
function plot_sizing_contour(Xpts, X_BE_spline)

    % compute the sizing function for each node on the mesh
    h = sizing(Xpts, X_BE_spline);

    %create contour plot
    x = Xpts(:, 1);
    y = Xpts(:, 2);
    [X, Y] = meshgrid(x,y);
    contour(X, Y, h);
    grid on
    title();

end


%Chat-GPT logic-checked and produced a cleaned version of this file
