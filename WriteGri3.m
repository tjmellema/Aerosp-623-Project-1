%% write_passage_coarse_gri_periodicShift.m
% Writes a .gri file from a PDE Toolbox mesh (msh) and automatically:
%  - extracts boundary edges from the mesh
%  - classifies boundary groups: PeriodicLeft, PeriodicRight, Bottom, Top
%  - writes nodes, boundary groups, element group
%  - builds periodic node pairs allowing a vertical shift (pitch-like shift)
%
% REQUIREMENTS BEFORE RUNNING:
%   msh : mesh object returned by generateMesh(model,...)
%   (optional) xL, xR : your geometry periodic plane locations (for printing only)
%
% OUTPUT:
%   passage_coarse.gri

clearvars -except msh xL xR; close all; clc;

assert(exist('msh','var')==1, 'Expected mesh "msh" in workspace (from generateMesh).');

%% ---------- BASIC MESH DATA ----------
nodeCoords = msh.Nodes;      % 2 x nNode
elemConn   = msh.Elements;   % 3 x nElemTot (linear triangles)

nNode    = size(nodeCoords, 2);
nElemTot = size(elemConn,   2);
Dim      = 2;

%% ---------- OPEN OUTPUT FILE ----------
outFile = 'passage_coarse.gri';
fid = fopen(outFile, 'w');
assert(fid > 0, 'Could not open output file: %s', outFile);

% Header: nNode nElemTot Dim
fprintf(fid, '%d %d %d\n', nNode, nElemTot, Dim);

% Node coordinates (space-separated)
for i = 1:nNode
    fprintf(fid, '%.16g %.16g\n', nodeCoords(1,i), nodeCoords(2,i));
end

%% ---------- EXTRACT BOUNDARY EDGES ----------
nodeXY     = nodeCoords.';   % nNode x 2
elemToNode = elemConn.';     % nElemTot x 3
TR = triangulation(elemToNode, nodeXY);
boundaryEdgesAll = freeBoundary(TR); % nBFaceTotal x 2

% Discover actual left/right boundary x-planes from mesh boundary nodes
boundaryNodeIds = unique(boundaryEdgesAll(:));
boundaryNodeX   = nodeXY(boundaryNodeIds, 1);

xLeftPlane  = min(boundaryNodeX);
xRightPlane = max(boundaryNodeX);

fprintf('Mesh boundary planes: xLeftPlane=%.16g, xRightPlane=%.16g\n', xLeftPlane, xRightPlane);
if exist('xL','var')==1 && exist('xR','var')==1
    fprintf('Your geometry planes: xL=%.16g, xR=%.16g\n', xL, xR);
end

% Classify periodic edges using midpoint-x (robust)
xSpan = xRightPlane - xLeftPlane;
xTol  = 1e-10 * max(1, abs(xSpan));     % if classification fails, loosen to 1e-9

edgeMidX = 0.5*( nodeXY(boundaryEdgesAll(:,1),1) + nodeXY(boundaryEdgesAll(:,2),1) );
isLeftEdge  = abs(edgeMidX - xLeftPlane)  < xTol;
isRightEdge = abs(edgeMidX - xRightPlane) < xTol;
isWallEdge  = ~(isLeftEdge | isRightEdge);

% Split wall edges into Bottom/Top using median of wall edge midpoint y
edgeMidY = 0.5*( nodeXY(boundaryEdgesAll(:,1),2) + nodeXY(boundaryEdgesAll(:,2),2) );
ySplit = median(edgeMidY(isWallEdge));

isBottomEdge = isWallEdge & (edgeMidY <  ySplit);
isTopEdge    = isWallEdge & (edgeMidY >= ySplit);

% Store boundary groups as unique undirected edges
bndGroups = cell(4,1);
bndGroups{1} = unique(sort(boundaryEdgesAll(isLeftEdge ,:),2), 'rows'); % PeriodicLeft
bndGroups{2} = unique(sort(boundaryEdgesAll(isRightEdge,:),2), 'rows'); % PeriodicRight
bndGroups{3} = unique(sort(boundaryEdgesAll(isBottomEdge,:),2),'rows'); % Bottom
bndGroups{4} = unique(sort(boundaryEdgesAll(isTopEdge   ,:),2),'rows'); % Top

bndTitles = {'PeriodicLeft','PeriodicRight','Bottom','Top'};
nBGroup = 4;

fprintf('Boundary faces total: %d\n', size(boundaryEdgesAll,1));
fprintf('Left=%d Right=%d Bottom=%d Top=%d\n', ...
    size(bndGroups{1},1), size(bndGroups{2},1), size(bndGroups{3},1), size(bndGroups{4},1));

%% ---------- WRITE BOUNDARY GROUPS ----------
fprintf(fid, '%d\n', nBGroup);

nf = 2; % edge has 2 nodes in 2D
for g = 1:nBGroup
    edges = bndGroups{g};
    fprintf(fid, '%d %d %s\n', size(edges,1), nf, bndTitles{g});
    for k = 1:size(edges,1)
        fprintf(fid, '%d %d\n', edges(k,1), edges(k,2));
    end
end

%% ---------- WRITE ELEMENT GROUP ----------
% One element group: all linear triangles
order = 1;
basis = 'TriLagrange';
fprintf(fid, '%d %d %s\n', nElemTot, order, basis);

for e = 1:nElemTot
    fprintf(fid, '%d %d %d\n', elemConn(1,e), elemConn(2,e), elemConn(3,e));
end

%% ---------- PERIODIC GROUP (WITH VERTICAL SHIFT) ----------
% Extract node sets on periodic boundaries
leftNodes  = unique(bndGroups{1}(:));
rightNodes = unique(bndGroups{2}(:));

fprintf('Periodic node counts: left=%d right=%d\n', numel(leftNodes), numel(rightNodes));
fprintf('Left y-range:  [%g, %g]\n', min(nodeXY(leftNodes,2)),  max(nodeXY(leftNodes,2)));
fprintf('Right y-range: [%g, %g]\n', min(nodeXY(rightNodes,2)), max(nodeXY(rightNodes,2)));

assert(numel(leftNodes) == numel(rightNodes), ...
    'Periodic node count mismatch: left=%d right=%d', numel(leftNodes), numel(rightNodes));

yLeft  = nodeXY(leftNodes,  2);
yRight = nodeXY(rightNodes, 2);

% Estimate vertical shift (robust): how much RIGHT must be shifted to match LEFT
yShift = median(yLeft) - median(yRight);
fprintf('Estimated vertical periodic shift yShift = %.6g\n', yShift);

% Pair by closest shifted-y: want yRight ≈ yLeft - yShift
periodicPairs = zeros(numel(leftNodes), 2);
usedRight = false(numel(rightNodes), 1);

% Tolerance for matching in y after shift (mesh nodes won't match to 1e-12)
yTol = 1e-3; % tighten later if you want

for a = 1:numel(leftNodes)
    targetY = yLeft(a) - yShift;
    [bestDist, idx] = min(abs(yRight - targetY));

    assert(bestDist < yTol, ...
        'No good shifted-y match for left node %d. Best |dy|=%g', leftNodes(a), bestDist);
    assert(~usedRight(idx), ...
        'Right node %d matched twice; periodic sets not one-to-one.', rightNodes(idx));

    usedRight(idx) = true;
    periodicPairs(a,:) = [leftNodes(a), rightNodes(idx)];
end

% Sort pairs by yLeft for nicer output
[~, ord] = sort(yLeft);
periodicPairs = periodicPairs(ord,:);

yMismatchShifted = max(abs( (nodeXY(periodicPairs(:,1),2) - yShift) - nodeXY(periodicPairs(:,2),2) ));
fprintf('Max |(yL - yShift) - yR|: %.3e\n', yMismatchShifted);

% Write periodic group header + contents
nPG = 1;
fprintf(fid, '%d %s\n', nPG, 'PeriodicGroup');

nPairs = size(periodicPairs, 1);
fprintf(fid, '%d %s\n', nPairs, 'Translational');

for k = 1:nPairs
    fprintf(fid, '%d %d\n', periodicPairs(k,1), periodicPairs(k,2));
end

%% ---------- CLOSE ----------
fclose(fid);
fprintf('Wrote %s successfully.\n', outFile);

%% ---------- OPTIONAL VISUAL CHECK ----------
figure('Name','Periodic nodes'); hold on; axis equal; grid on;
scatter(nodeXY(leftNodes,1),  nodeXY(leftNodes,2), 30, 'filled');
scatter(nodeXY(rightNodes,1), nodeXY(rightNodes,2), 30, 'filled');
legend('PeriodicLeft nodes','PeriodicRight nodes');
xlabel('x'); ylabel('y');
title(sprintf('Periodic node sets (estimated yShift = %.4g)', yShift));