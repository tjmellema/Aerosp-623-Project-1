%finds all boundary edges to avoid for smoothing
%
%B2E = matrix that stores the elemetns for each global boundary edge
%E2N = stores 3 nodes for each element
%num_nodes = number of nodes in the mesh
function isBoundary = BuildBoundaryNodes(B2E, E2N, num_nodes)

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
