import sizing_fucntion.m

%creates matrix that stores global edge number for every global element number
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

%TODO creates a set of node connectivities between each node
%
%E2N = element_data = matrix that stores 3 nodes for each element
function N = connectivity(E2N)
    % N{i} = list of neighboring nodes of node i
    N = cell(num_nodes,1);
    
    %for every element
    for e = 1:size(E2N,1)
    
        tri = E2N(e,:);
    
        n1 = tri(1);
        n2 = tri(2);
        n3 = tri(3);
    
        % Each node connects to the other two
        N{n1} = [N{n1}, n2, n3];
        N{n2} = [N{n2}, n1, n3];
        N{n3} = [N{n3}, n1, n2];
    end
    
    % Remove duplicates
    for i = 1:num_nodes
        N{i} = unique(N{i});
    end
    %TODO - Handle periodic boundaries
end %generated with Chat-GPT

%finds the midpoint location cordinate of an edge between two node cords
%
%Xpt1 = cordinate of node 1
%Xpt2 = cordinate of node 2
function X_edge_midpt = edge_midpoint(Xpt1, Xpt2)
    X_edge_midpt = (Xpt1 + Xpt2) ./ 2;
end

%takes a edge length and compares it to the sizing function
%
%node_data = point cordinates of each global node
%IE = nodes between global interior edges
%BE = nodes between global boundary edges
%X_BE_spline = splined boundary edge used for the shape function
function [Iflagged, Bflagged] = flag_edge(node_data, IE, BE, X_BE_spline)
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

    h_desired = sizing_function.sizing(X_edge_midpt, X_BE_spline);

    if h_desired < h_edge
        flagged = true;
    end
    %TODO
end

%average the node location with its adjacent nodes
%
%Xs = node_data = coordinate data for each node
%N = connectivity matrix
%isBoundary = boolean area for each node and if that node is a boundary
%num_iter = number of times to smooth the function
function Xsnew = smooth(Xs, N, isBoundary, num_iter)
    
    w = 0.3; % relaxation factor
    Xsnew = Xs; % initialize
    %loop through the smoothing function for each iteration
    for n = 1:num_iter
        %loop through every point
        for i = 1:size(Xs,1)
            %ensure we aren't smoothing boundary nodes
            if ~isBoundary(i)
        
                xi = Xs(i,:);
                avg = mean(Xs(N{i},:),1);
        
                Xsnew(i,:) = (1-w)*xi + w*avg;
        
            end
        end
        Xs = Xsnew;
    end
end

%finds all boundary edges to avoid for smoothing
%
%B2E = matrix that stores the elemetns for each global boundary edge
%E2N = stores 3 nodes for each element
%num_nodes = number of nodes in the mesh
function isBoundary = build_boundary_nodes(B2E, E2N, num_nodes)

    isBoundary = false(num_nodes,1);
    
    % local face → local node indices
    faceNodes = [2 3;
                 3 1;
                 1 2];
    
    for k = 1:size(B2E,1)
    
        elem = B2E(k,1);
        face = B2E(k,2);
    
        % global nodes of this element
        tri = E2N(elem,:);
    
        % nodes on this boundary edge
        n1 = tri(faceNodes(face,1));
        n2 = tri(faceNodes(face,2));
    
        isBoundary(n1) = true;
        isBoundary(n2) = true;
    end
end %written by Chat-GPT

function insertNode()
    %TODO
end

function createEdge()
    %TODO
end