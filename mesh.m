clc
clear all 

BladeLower = load('bladeupper.txt');  % lower surface per PDF
BladeUpper = load('bladelower.txt');  % Upper surface from PDF

figure; hold on; axis equal
plot(BladeLower(:,1), BladeLower(:,2), 'r.-')
plot(BladeUpper(:,1), BladeUpper(:,2), 'b.-')
legend('Lower','Upper')


plot(BladeLower(1,1), BladeLower(1,2), 'ko', 'MarkerFaceColor','k')
plot(BladeLower(end,1), BladeLower(end,2), 'ks', 'MarkerFaceColor','k')


% Make one closed polygon around the blade:
% Build closed polygon: surfaces are already oriented oppositely in the data
blade = [BladeUpper; BladeLower];

% Ensure it is explicitly closed (first point = last point)
if norm(blade(1,:) - blade(end,:)) > 1e-12
    blade = [blade; blade(1,:)];
end

figure; hold on; axis equal
plot(blade(:,1), blade(:,2), 'm.-');
title('Closed blade polygon (piecewise linear)');

% Find the xmin, xmax, ymin, ymax of the outside rectangle
xmin = min(blade(:, 1)); xmax = max(blade(:, 1));
ymin = min(blade(:, 2)); ymax = max(blade(:, 2));

% define a chord length to add padding to the rectangle
c = xmax - xmin;
% add padding
x1 = xmin - 0.5*c;
x2 = xmax + 1*c;
y1 = ymin - 0.5*c;
y2 = ymax + 0.5*c;


% Create a rectangle around the blade polygon
rectangle('Position', [x1, y1, x2 - x1, y2 - y1], 'EdgeColor', 'g', 'LineStyle', '--');
xlabel('X Coordinate');
ylabel('Y Coordinate');

%% BUILD PDE Geometry
% FIRST STEP Create a PDE model for the blade geometry
model = createpde();

% SECOND STEP Rectangle geometry matrix
R = [2, 4, x1, x2, x2, x1, y1, y1, y2, y2]';

% Add rectangle geometry to the PDE model
%geometryFromEdges(model, R1);

% THIRD STEP blade geometry matrix
xb = blade(:, 1);
yb = blade(:, 2);

% Remove any duplicate points (decsg often prefers open polygon list)
% 1) Remove explicit closure point if present (last == first)
if hypot(xb(end)-xb(1), yb(end)-yb(1)) < 1e-12
    xb(end) = [];
    yb(end) = [];
end

% 2) Remove consecutive duplicate points (common at the join)
d = hypot(diff(xb), diff(yb));
keep = [true; d > 1e-12];
xb = xb(keep);
yb = yb(keep);

% 3) Remove any remaining duplicates anywhere (rare but possible)
P = [xb yb];
[~, ia] = unique(P, 'rows', 'stable');
xb = xb(ia);
yb = yb(ia);

B = [2; numel(xb); xb(:); yb(:)]; % creates blade geometry matrix as a polygon sith xb sides

n = max(numel(R, B)); % find the max elements in the vector describing the geometry to pad them
% set the geometry description
geomDesc = [padcol(R, n), padcol(B, n)]; % in theory we don't know which one to pad so just pad both

% FOURTH STEP use set formula
ns = char('R1', 'B1')'; % naming convention from matlab tutorial
sf = 'R1 - B1'; % set formula

% FIFTH STEP build geometry using dl = decsg()
dl = decsg(geomDesc, sf, ns);

% FINAL STEP output pde geometry with edge labels
geometryFromEdges(model, dl);

% plot geometry
figure;
pdegplot(model,'EdgeLabels','on'); axis equal
title('PDE geometry with edge labels');


%% GENERATE COARSE MESH
Hmax = 0.1 * c; % just using a guess and tuning
msh = generateMesh(model,'Hmax',Hmax,'GeometricOrder','linear');

figure; pdemesh(model); axis equal
title(sprintf('Hmax=%.4g, Triangles=%d', Hmax, size(msh.Elements,2)));

% padcol function adds zeros to end of matrix to make them the same size
function c = padcol(v, n)
v = v(:);
c = [v; zeros(n-numel(v),1)];
end