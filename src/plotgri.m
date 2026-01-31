function plotgri(grifile);

% Read mesh (gri file)

fid = fopen(grifile, 'r');
% Read in nodes
A = fscanf(fid,'%d', 3);
nnode    = A(1);
nelemtot = A(2);
dim      = A(3);
V = zeros(nnode, dim);
for inode = 1:nnode,
  A = fscanf(fid, '%lf', 2);
  V(inode,:) = A(1:2)';
end
% Read through boundary info
figure(1); clf; hold on;
A = fscanf(fid, '%d', 1);
nbfgrp = A(1);
colors = lines(nbfgrp);
for ibfgrp = 1:nbfgrp
    A = textscan(fid,'%d %d %s', 1);
    nbface = A{1};
    nnode = A{2};
    title = A{3};
  for ibface = 1:nbface,
    A = fscanf(fid, '%d', nnode);
    plotedge_color(V(A,:), colors(ibfgrp,:));
  end
end

% Read in elements and plot edges
curtot = 0;
E2N = zeros(nelemtot, 6);
while (curtot ~= nelemtot),
  fgets(fid);
  sline = fgets(fid);
  [nelem, p, sbasis] = strread(sline, '%d %d %s');
  switch sbasis{1}
    case 'TriLagrange'
      nnode = (p+1)*(p+2)/2;
      nedge = 3;
      fvec = zeros(3,p+1);
      fvec(1,:) = [1:(p+1)];
      v = p+1; d = p;
      for k=1:(p+1), fvec(2,k) = v; v = v+d; d = d-1; end;
      v = 1; d = p+1;
      for k=1:(p+1), fvec(3,k) = v; v = v+d; d = d-1; end;
    otherwise
      error('element type not understood');
  end
  
  for elem = 1:nelem,
    A = fscanf(fid, '%d', nnode);  % HO nodes included
    for edge=1:nedge,
      plotedge(V(A(fvec(edge,:)),:));
    end
  end
  curtot = curtot + nelem;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified function such that periodic pairs area plotted with a red-green
% Also plots the node ids.
% dot for the pair
% Author: landinjm
% Begin edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
periodic_pairs = zeros(0, 2);
A = textscan(fid,'%d %s', 1);
n_periodic_groups = A{1};
assert(strcmp(A{2}, 'PeriodicGroup'), "Invalid .gri syntax for periodicity")
for i = 1:n_periodic_groups
    A = textscan(fid,'%d %s', 1);
    n_periodic_node_pairs = A{1};
    periodicity_type = A{2};
    assert(strcmp(periodicity_type,'Translational'), "Rotational periodicity not implemented")

    for j = 1:n_periodic_node_pairs
        periodic_pair = fscanf(fid, '%d', 2);
        periodic_pairs(end+1, :) = periodic_pair;
    end
end
scatter(V(periodic_pairs(:, 1),1), V(periodic_pairs(:, 1),2), 50, 'r', 'filled')
scatter(V(periodic_pairs(:, 2),1), V(periodic_pairs(:, 2),2), 50, 'g', 'filled')

for i = 1:size(V,1)
    text(V(i,1), V(i,2), sprintf('%d', i), 'FontSize', 8, 'Color', 'b', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fid);
nelem = nelemtot;
axis equal
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
box off; axis off
set(gca, 'LineWidth', 0.5);

%------------------------------------------------
function plotedge(V);
x = V(:,1);
y = V(:,2);
plot(x,y, 'k-', 'linewidth', 1);

function plotedge_color(V, color);
    x = V(:,1);
    y = V(:,2);
    plot(x, y, '-', 'Color', color, 'LineWidth', 4);


