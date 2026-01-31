
%plot the sizing function
[node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("passage_coarse.gri");

top_nodes = unique(boundary_mappings(boundary_mappings(:,3)==7,1:2));
bot_nodes = unique(boundary_mappings(boundary_mappings(:,3)==8,1:2));

top_coords = node_data(top_nodes,:);
bot_coords = node_data(bot_nodes,:);

X_BE_spline = spline_boundary(top_coords, bot_coords, 3);

plot_sizing_contour(node_data, X_BE_spline);

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
