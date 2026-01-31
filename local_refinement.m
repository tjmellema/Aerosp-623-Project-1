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

        edge_length = sqrt(sum((left_element_positions(1, :) - left_element_positions(2, :)).^2));

        % Grab the mesh length from the size function
        ideal_h = sizing(face_center, top_spline, bottom_spline);

        % Add the refinement to the left and right elements to make
        % things easier later on
        if edge_length > ideal_h
            edge_refinement(left_element, left_face) = 1;
            edge_refinement(right_element, right_face) = 1;
            interior_edge_refinement(i) = 1;
        end
    end

    % Initialize new nodes and elements
    edge2node = containers.Map('KeyType','char','ValueType','int32');
    edge_key = @(a,b) sprintf('%d_%d', min(a,b), max(a,b));
    new_nodes = nodes;
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
            eR = I2E(i,3);
            fR = I2E(i,4);

            nodesL = elements(eL,:);
            nodesL(fL) = [];

            nodesR = elements(eR,:);
            nodesR(fR) = [];

            a = nodesL(1);
            b = nodesL(2);

            c = nodesR(1);
            d = nodesR(2);
        
            key = edge_key(a,b);
            keyR = edge_key(c,d);
            if ~isKey(edge2node, key)
                midpoint_L = 0.5*(nodes(a,:) + nodes(b,:));
                midpoint_R = 0.5*(nodes(c,:) + nodes(d,:));
                
                if midpoint_L == midpoint_R
                    new_nodes(end+1,:) = midpoint_L;
                    edge2node(key) = size(new_nodes,1);
                else
                    new_nodes(end+1,:) = midpoint_L;
                    edge2node(key) = size(new_nodes,1);
                    new_nodes(end+1,:) = midpoint_R;
                    edge2node(keyR) = size(new_nodes,1);
                end
                
            end
        end
    end

    new_elements = zeros(0,3);

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

        flags = [~isempty(m12), ~isempty(m23), ~isempty(m31)];
    
        nflag = sum(flags);
    
        switch nflag
            % TODO: Double check these
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
                % Coordinates
                p1 = nodes(n1,:);
                p2 = nodes(n2,:);
                p3 = nodes(n3,:);
            
                % Interior angles
                a2 = triAngle(p1,p2,p3);
                a3 = triAngle(p2,p3,p1);
                a1 = triAngle(p3,p1,p2);
            
                % Find largest angle
                [~, idx] = max([a1 a2 a3]);
            
                if ~flags(1)   % edge 23 unrefined
                    % edges 12 and 31 refined → split across angle at node 2
                    disp('edge 23');
                    if idx == 2
                        % optimal split
                        new_elements(end+1,:) = [n1 m12 m31];
                        new_elements(end+1,:) = [m12 n2 m23];
                        new_elements(end+1,:) = [m31 m23 n3];
                    else
                        % default split
                        new_elements(end+1,:) = [n1 n2 m23];
                        new_elements(end+1,:) = [n1 m23 m31];
                        new_elements(end+1,:) = [m31 m23 n3];
                    end
            
                elseif ~flags(2)   % edge 31 unrefined
                    disp('edge 31');
                    % edges 12 and 23 refined → split across angle at node 3
                    if idx == 3
                        % optimal split
                        new_elements(end+1,:) = [n1 m12 m31];
                        new_elements(end+1,:) = [m12 n2 m23];
                        new_elements(end+1,:) = [m31 m23 n3];
                    else
                        % default split
                        new_elements(end+1,:) = [n1 m12 n3];
                        new_elements(end+1,:) = [m12 n2 m31];
                        new_elements(end+1,:) = [m31 n2 n3];
                    end
            
                else   % edge 12 unrefined
                    disp('edge 12');
                    % edges 23 and 31 refined → split across angle at node 1
                    if idx == 1
                        % optimal split
                        new_elements(end+1,:) = [n1 m12 m31];
                        new_elements(end+1,:) = [m12 n2 m23];
                        new_elements(end+1,:) = [m31 m23 n3];
                    else
                        % default split
                        new_elements(end+1,:) = [n1 m12 m23];
                        new_elements(end+1,:) = [m12 n2 m23];
                        new_elements(end+1,:) = [n1 m23 n3];
                    end
                end
            case 3
                new_elements(end+1,:) = [n1 m12 m31];
                new_elements(end+1,:) = [m12 n2 m23];
                new_elements(end+1,:) = [m31 m23 n3];
                new_elements(end+1,:) = [m12 m23 m31];
        end
    end
end

function ang = triAngle(A, B, C)
% TRIANGLE Returns the interior angle at vertex A (radians)
%
%   A, B, C are 1×2 coordinate vectors

    v1 = B - A;
    v2 = C - A;

    % Numerically safe acos
    cosang = dot(v1, v2) / (norm(v1) * norm(v2));
    cosang = max(min(cosang,1),-1);   % clamp for safety

    ang = acos(cosang);
end
