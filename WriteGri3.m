%% write_passage_coarse_gri.m
% Writes a .gri file from the nodeCoords and elemConn arrays fixed in meshWithoutGri
% classification: PeriodicLeft, PeriodicRight, Wall (Bottom/Top)

% FIX 1: If msh was cleared but nodeCoords/elemConn exist, use them.
if ~exist('msh','var')
    if exist('nodeCoords','var') && exist('elemConn','var')
        fprintf('msh not found, using nodeCoords and elemConn from workspace.\n');
    else
        error('Required mesh data (msh or nodeCoords/elemConn) not found.');
    end
else
    nodeCoords = msh.Nodes;
    elemConn   = msh.Elements;
end

nNode    = size(nodeCoords, 2);
nElemTot = size(elemConn,   2);
nodeXY   = nodeCoords.'; 
Dim      = 2;
pitch    = 18; % Required for pairing

%% ---------- OPEN OUTPUT FILE ----------
outFile = 'passage_coarse.gri';
fid = fopen(outFile, 'w');
fprintf(fid, '%d %d %d\n', nNode, nElemTot, Dim);
for i = 1:nNode
    fprintf(fid, '%.16g %.16g\n', nodeCoords(1,i), nodeCoords(2,i));
end

%% ---------- EXTRACT BOUNDARY EDGES ----------
elemToNode = elemConn.';
TR = triangulation(elemToNode, nodeXY);
boundaryEdgesAll = freeBoundary(TR); 

edgeMidX = 0.5*( nodeXY(boundaryEdgesAll(:,1),1) + nodeXY(boundaryEdgesAll(:,2),1) );
edgeMidY = 0.5*( nodeXY(boundaryEdgesAll(:,1),2) + nodeXY(boundaryEdgesAll(:,2),2) );

xLeftPlane  = min(nodeXY(:,1));
xRightPlane = max(nodeXY(:,1));
yBottom     = min(nodeXY(:,2));
yTop        = max(nodeXY(:,2));

% Classify boundary groups
isLeftEdge   = abs(edgeMidX - xLeftPlane)  < 1e-8;
isRightEdge  = abs(edgeMidX - xRightPlane) < 1e-8;
isBottomEdge = abs(edgeMidY - yBottom)     < 1e-8;
isTopEdge    = abs(edgeMidY - yTop)        < 1e-8;

bndGroups = {
    unique(sort(boundaryEdgesAll(isLeftEdge,:),2), 'rows');  % 1: Inlet
    unique(sort(boundaryEdgesAll(isRightEdge,:),2), 'rows'); % 2: Outlet
    unique(sort(boundaryEdgesAll(isBottomEdge,:),2), 'rows');% 3: Bottom (Periodic)
    unique(sort(boundaryEdgesAll(isTopEdge,:),2), 'rows');   % 4: Top (Periodic)
};
bndTitles = {'Inlet','Outlet','Bottom','Top'};

fprintf(fid, '%d\n', 4);
for g = 1:4
    edges = bndGroups{g};
    fprintf(fid, '%d 2 %s\n', size(edges,1), bndTitles{g});
    for k = 1:size(edges,1)
        fprintf(fid, '%d %d\n', edges(k,1), edges(k,2));
    end
end

%% ---------- WRITE ELEMENTS ----------
fprintf(fid, '%d 1 TriLagrange\n', nElemTot);
for e = 1:nElemTot
    fprintf(fid, '%d %d %d\n', elemConn(1,e), elemConn(2,e), elemConn(3,e));
end

%% ---------- PERIODIC GROUP (TOP-BOTTOM) ----------
% According to Section 1.2 of proj.pdf, nodes on periodic boundaries match.
% In this blade passage, the Bottom wall and Top wall are the periodic pairs.
botNodes = unique(bndGroups{3}(:));
topNodes = unique(bndGroups{4}(:));

periodicPairs = [];
for i = 1:length(botNodes)
    cX = nodeXY(botNodes(i), 1);
    cY = nodeXY(botNodes(i), 2);
    
    % Match node on top boundary with same X and Y + pitch
    dist = sqrt((nodeXY(topNodes,1) - cX).^2 + (nodeXY(topNodes,2) - (cY + pitch)).^2);
    [minD, idx] = min(dist);
    
    if minD < 1e-4
        periodicPairs = [periodicPairs; botNodes(i), topNodes(idx)];
    end
end

fprintf(fid, '1 PeriodicGroup\n');
fprintf(fid, '%d Translational\n', size(periodicPairs, 1));
for k = 1:size(periodicPairs, 1)
    fprintf(fid, '%d %d\n', periodicPairs(k,1), periodicPairs(k,2));
end

fclose(fid);
fprintf('Successfully wrote %s\n', outFile);