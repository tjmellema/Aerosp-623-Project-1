%takes a edge length and compares it to the sizing function
%
%node_data = point cordinates of each global node
%IE = nodes between global interior edges
%BE = nodes between global boundary edges
%X_BE_spline = splined boundary edge used for the shape function
function [Iflagged, Bflagged] = FlagEdge(node_data, IE, BE, X_BE_spline)
    %initialize the flagged array for every edge
    Iflagged = false(size(IE, 1));
    Bflagged = false(size(BE, 1));

    %loop over the interior edges
    for i = 1:length(IE)
        %get the nodes that make up an edge 
        [n1, n2, elem1, elem2] = IE(i);
        Xpt1 = node_data(n1); Xpt2 = node_data(n2);

        %get the lengths of a current edge
        h_current = norm(Xpt2 - Xpt1);
        X_midpt = edge_midpoint(Xpt1, Xpt2);
        h_desired = sizing_function.sizing(X_midpt, X_BE_spline);

        %ensure the edge length is below the needed amount
        if h_desired < h_current
            Iflagged(i) = true;
        end
    end

    %loop over the boundary edges
    for i = 1:length(BE)
        %get the nodes that make up an edge 
        [n1, n2, elem] = BE(i);
        Xpt1 = node_data(n1); Xpt2 = node_data(n2);

        %get the lengths of a current edge
        h_current = norm(Xpt2 - Xpt1);
        X_midpt = edge_midpoint(Xpt1, Xpt2);
        h_desired = sizing_function.sizing(X_midpt, X_BE_spline);

        %ensure the edge length is below the needed amount
        if h_desired < h_current
            Bflagged(i) = true;
        end

    end
end

%finds the midpoint location cordinate of an edge between two node cords
%
%Xpt1 = cordinate of node 1
%Xpt2 = cordinate of node 2
function X_edge_midpt = edge_midpoint(Xpt1, Xpt2)
    X_edge_midpt = (Xpt1 + Xpt2) ./ 2;
end