import sizing_fucntion.m
import mesh.m


%creates matrix that stores global edge number for every global node number
%
%I2E = interior edge matrix storing each local edge number for each node
%B2E = boundary edge matrix storing each local edge number for each node
%num_nodes =
function E2I = compute_E2I(I2E, B2E, num_nodes)
    E2I = zeros(num_nodes, 3);
    %loop over I2E
    for i = 1:length(I2E)
        [elemL, faceL, elemR, faceR] = I2E(i);
        E2I(elemL, faceL) = i;
        E2I(elemR, faceR) = i;
    end

    %loop over B2E
    for i = 1:length(BE)
        [elem, face, bgroup] = B2E(i);
        E2I(elem, face) = -1 .* i;
    end
end

%finds the midpoint location cordinate of an edge between two node cords
%
%Xpt1 = cordinate of node 1
%Xpt2 = cordinate of node 2
function X_edge_midpt = edge_midpoint(Xpt1, Xpt2)
    X_edge_midpt = (Xpt1 + Xpt2) ./ 2;
end

%takes a edge length and compares it to the sizing function
%
%X_edge_midpt = midpoint coordinate of hte desired edge
%h_edge = current edge length of desired edge
%X_BE_spline = splined boundary edge used for the shape function
function flagged = flag_edge(X_edge_midpt, h_edge, X_BE_spline)
    flagged = false;
    h_desired = sizing_function.sizing(X_edge_midpt, X_BE_spline);

    if h_desired < h_edge
        flagged = true;
    end
    %TODO
end

%average the node location with its adjacent nodes
%
%X = current node to average out
function Xnew = smooth(X)
    w = 0.3; %weight factor
    Xnew = (1-w) .* X + w .* avg;
    %TODO
end