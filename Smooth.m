%average the node location with its adjacent nodes
%
%Xs = node_data = coordinate data for each node
%N = connectivity matrix
%isBoundary = boolean area for each node and if that node is a boundary
%num_iter = number of times to smooth the function
function Xsnew = Smooth(Xs, num_iter, E2N, num_nodes)
    %call the helper funcitons
    N = Connectivity(E2N, num_nodes); % Create connectivity matrix
    isBoundary = check_boundary(E2N, num_nodes); % Identify boundary nodes
    %{
    figure;
    scatter(Xs(:,1), Xs(:,2), 20, isBoundary, 'filled');
    axis equal; colorbar;
    title('Detected Boundary Nodes');
    %}
    
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
function isBoundary = check_boundary(E2N, num_nodes)

    edges = [E2N(:,[1 2]);
             E2N(:,[2 3]);
             E2N(:,[3 1])];

    edges = sort(edges,2);

    [u,~,idx] = unique(edges,'rows');
    counts = accumarray(idx,1);

    boundary_edges = u(counts == 1,:);

    isBoundary = false(num_nodes,1);
    isBoundary(boundary_edges(:)) = true;
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