function mesh_verification()
    % checks the closed-loop normal summation property (Geometric Conservation)
    clear; clc;
    
    % 1. Run Test on 'test.gri'
    fprintf('------------------------------------------------\n');
    fprintf('Running Verification on "test.gri"...\n');
    verify_mesh('test.gri');
    
    % % 2. Run Test on 'task1_mesh.gri' (whatever we call mesh.gri)
    % fprintf('\n------------------------------------------------\n');
    % fprintf('Running Verification on "task1_mesh.gri"...\n');
    % verify_mesh('task1_mesh.gri');
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
    % This includes periodic boundaries treated as interior edges
    for k = 1:size(I2E, 1)
        eL = I2E(k, 1); % Left element index
        eR = I2E(k, 3); % Right element index
        
        % In(k,:) contains the normal scaled by length (nx*l, ny*l)
        % Normal points from L to R
        vec_normal_l = In(k, :); 
        
        % Add to Left: vec_normal points OUT of L 
        Elem_Sum(eL, :) = Elem_Sum(eL, :) + vec_normal_l;
        
        % Subtract from Right: vec_normal points INTO R, so negative points OUT 
        Elem_Sum(eR, :) = Elem_Sum(eR, :) - vec_normal_l;
    end
    
    % Process Boundary Faces (Bn)
    for k = 1:size(B2E, 1)
        eL = B2E(k, 1); % Element adjacent to boundary
        
        % Bn(k,:) contains the outward normal scaled by length
        vec_normal_l = Bn(k, :);
        
        % Add to element: normal already points OUT of the domain
        Elem_Sum(eL, :) = Elem_Sum(eL, :) + vec_normal_l;
    end
    
    % 4. Compute Magnitude of Error Ee for each element
    Ee_mag = sqrt(Elem_Sum(:,1).^2 + Elem_Sum(:,2).^2);
    
    % 5. Results
    max_err = max(Ee_mag);
    avg_err = mean(Ee_mag);
    
    fprintf('  Number of Elements: %d\n', nElem);
    
    % Check for machine precision 
    if max_err < 1e-12
        fprintf('SUCCESS. Geometric Conservation Law satisfied.\n');
    else
        fprintf('FAIL. Check normal directions or periodic pairings.\n');
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