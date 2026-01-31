function isBoundary = BuildBoundaryNodes(BE, NE, num_nodes)

isBoundary = false(num_nodes,1);

for i = 1:size(BE,1)
    isBoundary(BE(i,1)) = true;
    isBoundary(BE(i,2)) = true;
end

end