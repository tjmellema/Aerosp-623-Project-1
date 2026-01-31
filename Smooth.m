%average the node location with its adjacent nodes
%
%Xs = node_data = coordinate data for each node
%N = connectivity matrix
%isBoundary = boolean area for each node and if that node is a boundary
%num_iter = number of times to smooth the function
function Xsnew = Smooth(Xs, num_iter, E2N, B2E, num_nodes)
    %call the helper funcitons
    N = Connectivity(E2N, num_nodes); % Create connectivity matrix
    isBoundary = check_boundary(E2N, B2E, num_nodes); % Identify boundary nodes


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

%create a of each node that is considered a boundary node
function isBoundary = check_boundary(E2N, B2E, num_nodes)
    isBoundary = false(num_nodes,1);
    for i = 1:size(B2E, 1)
        [elem, edge, ~] = B2E(i);

        %if this is local edge 1 - mark nodes 2 and 3
        if edge == 1
            isBoundary(E2N(elem, 2)) = true;
            isBoundary(E2N(elem, 3)) = true;
        
       %if this is local edge 2 - mark nodes 1 and 3
        elseif edge == 2
            isBoundary(E2N(elem, 1)) = true;
            isBoundary(E2N(elem, 3)) = true;

        %if this is local edge 3 - mark nodes 1 and 2
        else
            isBoundary(E2N(elem, 1)) = true;
            isBoundary(E2N(elem, 2)) = true;
        end
    end
end

%creates a set of node connectivities between each node
%
%E2N = element_data = matrix that stores 3 nodes for each element
function N = Connectivity(E2N, num_nodes)
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
end %generated with Chat-GPT