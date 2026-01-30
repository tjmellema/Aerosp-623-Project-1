function mesh_verification(filename)
    fprintf('------------------------------------------------\n');
    fprintf('Running Verification on "%s"...\n', filename);
    verify_mesh(filename);
end

function verify_mesh(filename)
    % 1. Read Mesh using your minimal reader
    [nodes, elem, b_groups, periodic_pairs] = read_gri(filename);
    nElem = size(elem, 1);
    
    % 2. Call the wrapper in 'mappings.m' to get outputs
    [I2E, B2E, In, Bn, Area] = mappings(nodes, elem, b_groups, periodic_pairs);
    
    % 3. Initialize Error Sum vectors for every element
    % E_e should be a vector quantity (nx*l, ny*l) 
    Elem_Sum = zeros(nElem, 2);
    
    % Process Interior Faces (In)
    for k = 1:size(I2E, 1)
        eL = I2E(k, 1); 
        eR = I2E(k, 3);
        
        % Identify face nodes to calculate length (L)
        % Local face i is opposite local node i
        fL = I2E(k, 2);
        face_node_indices = [1, 2, 3];
        face_node_indices(fL) = []; % Remove the opposite node
        
        n1 = elem(eL, face_node_indices(1));
        n2 = elem(eL, face_node_indices(2));
        
        % Calculate edge length L
        L = sqrt((nodes(n2,1) - nodes(n1,1))^2 + (nodes(n2,2) - nodes(n1,2))^2);

        % SCALE the partner's unit normal by the length
        vec_normal_l = In(k, :) * L;
        
        Elem_Sum(eL, :) = Elem_Sum(eL, :) + vec_normal_l;
        Elem_Sum(eR, :) = Elem_Sum(eR, :) - vec_normal_l;
    end
    
    % Process Boundary Faces (Bn)
    for k = 1:size(B2E, 1)
        eL = B2E(k, 1);
        fL = B2E(k, 2);

        % Identify face nodes to calculate length (L)
        face_node_indices = [1, 2, 3];
        face_node_indices(fL) = [];

        n1 = elem(eL, face_node_indices(1));
        n2 = elem(eL, face_node_indices(2));

        % Calculate edge length L
        L = sqrt((nodes(n2,1) - nodes(n1,1))^2 + (nodes(n2,2) - nodes(n1,2))^2);

        % SCALE the unit normal (Bn) by the length
        vec_normal_l = Bn(k, :) * L;

        Elem_Sum(eL, :) = Elem_Sum(eL, :) + vec_normal_l;
    end
    
    % 4. Compute Magnitude of Error Ee for each element
    Ee_mag = sqrt(Elem_Sum(:,1).^2 + Elem_Sum(:,2).^2);
    
    % 5. Results
    max_err = max(Ee_mag);
    
    fprintf('  Number of Elements: %d\n', nElem);
    fprintf('  Max. Error: %d\n', max_err);
    
    % Check for machine precision 
    if max_err < 1e-12
        fprintf('SUCCESS. Geometric Conservation Law satisfied.\n');
    else
        fprintf('FAIL. Max Error = %.4f\n', max_err);
    end
end

% -------------------------------------------------------------------------
% Helper Functions (Preprocessing)
% -------------------------------------------------------------------------

function [nodes, elem, b_groups] = read_gri_minimal(filename)
    % A minimal GRI reader to support the task
    fid = fopen(filename, 'r');
    if fid == -1, error(['Cannot open ' filename]); end
    
    % Header
    vals = fscanf(fid, '%d', 3);
    nNode = vals(1); nElem = vals(2);
    nodes = fscanf(fid, '%f', [2, nNode])';
    
    % B Groups
    nBGroup = fscanf(fid, '%d', 1);
    b_groups = struct('faces', {});
    
    for i = 1:nBGroup
        % Read group header (Faces, Order, Name)
        % Use fgets to handle string properly
        line = fgetl(fid); while isempty(strtrim(line)), line=fgetl(fid); end
        % Handle cases where the line might be just after the count
        if sscanf(line, '%d', 1) == nBGroup, line = fgetl(fid); end 
        
        c = textscan(line, '%d %d %s');
        nBF = c{1};
        % Read faces
        faces = fscanf(fid, '%d', [2, nBF])';
        b_groups(i).faces = faces;
    end
    
    % Elements
    line = fgetl(fid); while isempty(strtrim(line)), line=fgetl(fid); end
    % Skip group header info
    elem = fscanf(fid, '%d', [3, nElem])';
    
    fclose(fid);
end
