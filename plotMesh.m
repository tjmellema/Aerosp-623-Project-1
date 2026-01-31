function plotMesh(E2N, nodes, elementType, showNodeNumbers)
% plotMesh - Plots a 2D finite element mesh
%
% Syntax: plotMesh(E2N, nodes, elementType, showNodeNumbers)
%
% Inputs:
%   E2N - Element-to-node connectivity matrix (elements x nodes per element)
%   nodes - Node coordinates (numNodes x 2) [x, y]
%   elementType - 'quad' for quadrilaterals, 'tri' for triangles (default: 'quad')
%   showNodeNumbers - true/false to show node numbers (default: false)

    if nargin < 3
        elementType = 'quad';
    end
    if nargin < 4
        showNodeNumbers = false;
    end

    figure;
    hold on;
    axis equal;
    grid on;
    xlabel('X'); ylabel('Y');
    title('Mesh Plot');

    numElems = size(E2N, 1);

    for e = 1:numElems
        elemNodes = E2N(e, :);
        coords = nodes(elemNodes, :);
        
        % Close the loop
        if strcmpi(elementType, 'quad')
            coords = [coords; coords(1,:)];
        elseif strcmpi(elementType, 'tri')
            coords = [coords; coords(1,:)];
        else
            error('elementType must be ''quad'' or ''tri''');
        end
        
        plot(coords(:,1), coords(:,2), 'b-', 'LineWidth', 1.5);
    end

    % Optionally plot node numbers
    if showNodeNumbers
        numNodes = size(nodes,1);
        for n = 1:numNodes
            text(nodes(n,1), nodes(n,2), num2str(n), 'Color', 'r', 'FontSize', 8, 'HorizontalAlignment', 'center');
        end
    end

    hold off;
end