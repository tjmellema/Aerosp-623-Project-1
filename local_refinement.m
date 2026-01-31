function [new_nodes, new_elements] = local_refinement(nodes, elements,...
                                       top_spline, bottom_spline, ...
                                       I2E, B2E)
    n_elements = size(elements, 1);
    n_nodes = size(nodes, 1);
    n_interior_faces = size(I2E, 1);

    edge_refinement = zeros(n_elements, 3);
    interior_edge_refinement = zeros(n_interior_faces, 1);
    % Loop over the interior faces to determine whether to flag refinement
    % or not
    for i = 1:n_interior_faces
        % Grab the elements
        left_element = I2E(i, 1);
        left_face = I2E(i, 2);
        right_element = I2E(i, 3);
        right_face = I2E(i, 4);

        % Grab edge length from the left element
        left_elements_nodes = elements(left_element, :);
        left_elements_nodes(left_face) = [];
        left_element_positions = nodes(left_elements_nodes, :);

        face_center = mean(left_element_positions(1, :) + left_element_positions(2, :), 1);

        edge_length = sqrt(sum(left_element_positions(1, :) - left_element_positions(2, :)).^2);

        % Grab the mesh length from the size function
        ideal_h = sizing(face_center, top_spline, bottom_spline);

        % Add the refinement to the left and right elements to make
        % things easier later on
        if edge_length < ideal_h
            edge_refinement(left_element, left_face) = 1;
            edge_refinement(right_element, right_face) = 1;
            interior_edge_refinement(i) = 1;
        end
    end

    % Initialize new nodes and elements
    edge2node = containers.Map('KeyType','char','ValueType','int32');
    new_nodes = nodes;
    edge_key = @(a,b) sprintf('%d_%d', min(a,b), max(a,b));
    function m = get_mid(a,b)
        key = sprintf('%d_%d', min(a,b), max(a,b));
        if isKey(edge2node,key)
            m = edge2node(key);
        else
            m = [];
        end
    end

    % Update the nodes and elements with the refinement
    for i = 1:n_interior_faces
        if interior_edge_refinement(i) == 1
            eL = I2E(i,1);
            fL = I2E(i,2);

            nodesL = elements(eL,:);
            nodesL(fL) = [];
    
            a = nodesL(1);
            b = nodesL(2);
    
            key = edge_key(a,b);
    
            if ~isKey(edge2node, key)
                midpoint = 0.5*(nodes(a,:) + nodes(b,:));
                new_nodes(end+1,:) = midpoint;
                edge2node(key) = size(new_nodes,1);
            end
        end
    end

    new_elements = [];

    for e = 1:n_elements
        nodes_e = elements(e,:);
        flags   = edge_refinement(e,:);
    
        n1 = nodes_e(1);
        n2 = nodes_e(2);
        n3 = nodes_e(3);
    
        % Look up midpoints if they exist
        m12 = get_mid(n1,n2);
        m23 = get_mid(n2,n3);
        m31 = get_mid(n3,n1);
    
        nflag = sum(flags);
    
        switch nflag
    
            case 0
                new_elements(end+1,:) = nodes_e;
    
            case 1
                if flags(1)
                    new_elements(end+1,:) = [n1 m12 n3];
                    new_elements(end+1,:) = [m12 n2 n3];
                elseif flags(2)
                    new_elements(end+1,:) = [n1 n2 m23];
                    new_elements(end+1,:) = [n1 m23 n3];
                else
                    new_elements(end+1,:) = [n1 n2 m31];
                    new_elements(end+1,:) = [m31 n2 n3];
                end
    
            case 2
                if ~flags(1)
                    new_elements(end+1,:) = [n1 n2 m23];
                    new_elements(end+1,:) = [n1 m23 m31];
                    new_elements(end+1,:) = [m31 m23 n3];
                elseif ~flags(2)
                    new_elements(end+1,:) = [n1 m12 n3];
                    new_elements(end+1,:) = [m12 n2 m31];
                    new_elements(end+1,:) = [m31 n2 n3];
                else
                    new_elements(end+1,:) = [n1 m12 m23];
                    new_elements(end+1,:) = [m12 n2 n3];
                    new_elements(end+1,:) = [m23 m31 n3];
                end
    
            case 3
                new_elements(end+1,:) = [n1 m12 m31];
                new_elements(end+1,:) = [m12 n2 m23];
                new_elements(end+1,:) = [m31 m23 n3];
                new_elements(end+1,:) = [m12 m23 m31];
        end
    end
end
