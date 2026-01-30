%% mesh2_fixed_fullscript.m
% Turbine blade passage mesh (Figure 1 style) using MATLAB PDE Toolbox
% KEY: Keep your original file mapping, but SWAP which curve is used as
%      bottom vs top wall:
%   - Bottom wall = Lower (UNSHIFTED)
%   - Top wall    = Upper shifted up by pitch
%
% Left/right boundaries are periodic planes placed Lper upstream/downstream.

clear; close all; clc;

%% ---- USER SETTINGS ----
pitch = 18;      % mm (vertical spacing)
Lper  = 17;      % mm (streamwise periodic padding upstream/downstream)
Hmax_scale = 1/15;  % adjust to get ~1000–1500 triangles

% IMPORTANT: keep your current mapping (per your note)
UpperFile = 'contrib/bladelower.txt';   % "Upper" dataset variable (might need to swap)
LowerFile = 'contrib/bladeupper.txt';   % "Lower" dataset variable (might need to swap)

tol = 1e-10;

%% ---- LOAD BLADE DATA ----
assert(isfile(UpperFile), "Can't find %s", UpperFile);
assert(isfile(LowerFile), "Can't find %s", LowerFile);

Upper = load(UpperFile);  % [x y]
Lower = load(LowerFile);  % [x y]
assert(size(Upper,2)==2 && size(Lower,2)==2, 'Blade files must be Nx2 arrays');

% Sort by x for interpolation/clipping
Upper = sortrows(Upper, 1);
Lower = sortrows(Lower, 1);

%% ---- DEFINE BLADE x-EXTENT + PERIODIC PLANES ----
xminB = min([Upper(:,1); Lower(:,1)]);
xmaxB = max([Upper(:,1); Lower(:,1)]);
chord = xmaxB - xminB;

xL = xminB - Lper;  % inlet periodic plane
xR = xmaxB + Lper;  % outlet periodic plane

%% ---- SWAP USAGE (this is the fix) ----
% Bottom wall = Lower (unshifted)
Bottom0 = Lower;

% Top wall = Upper shifted up by pitch
Top0 = Upper;
Top0(:,2) = Top0(:,2) + pitch;

%% ---- CLIP WALL CURVES TO [xminB, xmaxB] ---- probably not necessary
BottomWall = Bottom0(Bottom0(:,1) >= xminB & Bottom0(:,1) <= xmaxB, :);
TopWall0   = Top0(Top0(:,1)     >= xminB & Top0(:,1)     <= xmaxB, :);

BottomWall = sanitize_polyline(BottomWall, tol);
TopWall0   = sanitize_polyline(TopWall0, tol);

% For CCW boundary assembly, top should run right->left
TopWall = flipud(TopWall0);

%% ---- ENDPOINT y VALUES at xminB/xmaxB ----
yB_xmin = interp1(BottomWall(:,1), BottomWall(:,2), xminB, 'spline','extrap');
yB_xmax = interp1(BottomWall(:,1), BottomWall(:,2), xmaxB, 'spline','extrap');
yT_xmin = interp1(TopWall0(:,1),   TopWall0(:,2),   xminB, 'spline','extrap');
yT_xmax = interp1(TopWall0(:,1),   TopWall0(:,2),   xmaxB, 'spline','extrap');

% Extend straight to periodic planes at same y as blade endpoints
yB_L = yB_xmin;  yB_R = yB_xmax;
yT_L = yT_xmin;  yT_R = yT_xmax;

%% ---- ASSEMBLE PASSAGE BOUNDARY LOOP (CCW, should not self-intersect) ----
passage = [
    xL    yB_L;
    xminB yB_xmin;

    BottomWall;      % bottom wall left->right

    xmaxB yB_xmax;
    xR    yB_R;

    xR    yT_R;
    xmaxB yT_xmax;

    TopWall;         % top wall right->left

    xminB yT_xmin;
    xL    yT_L;

    xL    yB_L
];

passage = sanitize_polygon_for_decsg(passage, tol);

%% ---- PLOT PASSAGE ----
figure('Name','Passage boundary'); hold on; axis equal; grid on;
plot(passage(:,1), passage(:,2), 'b.-');
xlabel('x (mm)'); ylabel('y (mm)');
title('Passage boundary loop (after swap fix)');

%% ---- CREATE PDE GEOMETRY ----
model = createpde();
xb = passage(:,1);
yb = passage(:,2);
G = [2; numel(xb); xb(:); yb(:)];

dl = decsg(G);
geometryFromEdges(model, dl);

figure('Name','Edge labels');
pdegplot(model, 'EdgeLabels','on'); axis equal; grid on;
title('PDE geometry with edge labels');

%% ---- GENERATE COARSE MESH ----
Hmax = chord * Hmax_scale;
msh = generateMesh(model, 'Hmax', Hmax, 'GeometricOrder','linear');

figure('Name','Coarse mesh');
pdemesh(model); axis equal; grid on;
title(sprintf('Coarse mesh: Hmax=%.4g, Triangles=%d', Hmax, size(msh.Elements,2)));

fprintf('Hmax = %.6g\n', Hmax);
fprintf('Triangles = %d\n', size(msh.Elements,2));

%% ---- QUICK CHECKS ----
fprintf('Upstream length:   %.6g (target %g)\n', xminB - xL, Lper);
fprintf('Downstream length: %.6g (target %g)\n', xR - xmaxB, Lper);
fprintf('Pitch check at xminB (top-bottom): %.6g (target %g)\n', yT_xmin - yB_xmin, pitch);
fprintf('Pitch check at xmaxB (top-bottom): %.6g (target %g)\n', yT_xmax - yB_xmax, pitch);

%% ================= Helper functions =================
function P = sanitize_polyline(P, tol)
% Remove NaNs, consecutive duplicates, and duplicates anywhere (stable)
P = P(all(isfinite(P),2), :);

d = hypot(diff(P(:,1)), diff(P(:,2)));
P = P([true; d > tol], :);

Q = round(P ./ tol) * tol;
[~, ia] = unique(Q, 'rows', 'stable');
P = P(ia, :);

d = hypot(diff(P(:,1)), diff(P(:,2)));
P = P([true; d > tol], :);
end

function P = sanitize_polygon_for_decsg(P, tol)
% Remove repeats and ensure polygon is NOT explicitly closed.
P = P(all(isfinite(P),2), :);

if hypot(P(end,1)-P(1,1), P(end,2)-P(1,2)) < tol
    P(end,:) = [];
end

d = hypot(diff(P(:,1)), diff(P(:,2)));
P = P([true; d > tol], :);

Q = round(P ./ tol) * tol;
[~, ia] = unique(Q, 'rows', 'stable');
P = P(ia, :);

d = hypot(diff(P(:,1)), diff(P(:,2)));
P = P([true; d > tol], :);
end


