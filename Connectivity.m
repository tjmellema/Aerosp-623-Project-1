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
    %TODO - Handle periodic boundaries
end %generated with Chat-GPT
