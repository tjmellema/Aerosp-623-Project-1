%requested function in 1.4 to show the projection function
%
%Xpts = set of coordinates to project
%X_BE = unsplined boundary coordinates
function plot_projection()
    % variables used for fig
    Xpts = [-.05,.1;.66,.5;.25,.03125;1.5,.25;1,.33];
    X_BE = [0,0;.5,.125;1,.125;1.5,.0625;2,0];
    ref = 3;

    %boundary spline
    [X_BE_spline, ~] = spline_boundary(X_BE, X_BE, ref);

    %projection
    for i = 1:size(Xpts,1)
        [d(i), xb(i), yb(i)] = projection(Xpts(i,:), X_BE_spline);
    end

    %plot the reference spline
    plot(X_BE_spline(:,1),X_BE_spline(:,2), Linewidth=1);
    hold on

    %plot the points
    scatter(Xpts(:,1),Xpts(:,2));

    %plot the projection to the spline
    X_plot = [Xpts(:,1)';xb];
    Y_plot = [Xpts(:,2)';yb];
    plot(X_plot, Y_plot, '-b', Linewidth=2);

    title('Sample Projection Function')
    grid on
    axis equal
    xlabel('x')
    ylabel('y')


end