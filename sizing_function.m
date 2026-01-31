
%plot the sizing function
[node_data, element_data, boundary_mappings, periodic_pairs] = read_gri("passage_coarse.gri");

top_nodes = unique(boundary_mappings(boundary_mappings(:,3)==7,1:2));
bot_nodes = unique(boundary_mappings(boundary_mappings(:,3)==8,1:2));

top_coords = node_data(top_nodes,:);
bot_coords = node_data(bot_nodes,:);

[X_BE_spline_top, X_BE_spline_bottom] = spline_boundary(top_coords, bot_coords, 3);

plot_sizing_contour(node_data, X_BE_spline_top, X_BE_spline_bottom);

%creates a contour of the plot
%
%Xpts = node coordinates of a mesh
%X_BE_spline = coordinates of the splined boundary used for sizing function
function plot_sizing_contour(Xpts, X_BE_spline_top, X_BE_spline_bottom)

    % compute the sizing function for each node on the mesh
    for i = 1:size(Xpts,1)
        h(i) = sizing(Xpts(i,:), X_BE_spline_top, X_BE_spline_bottom);
    end

    %create contour plot
    x = Xpts(:, 1);
    y = Xpts(:, 2);
    [X, Y] = meshgrid(x,y);
    Z = griddata(x, y, h', X, Y, 'natural');
    surf(x, y, h');
    grid on
    title();

end


%Chat-GPT logic-checked and produced a cleaned version of this file
