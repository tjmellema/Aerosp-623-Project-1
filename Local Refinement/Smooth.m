%average the node location with its adjacent nodes
%
%Xs = node_data = coordinate data for each node
%N = connectivity matrix
%isBoundary = boolean area for each node and if that node is a boundary
%num_iter = number of times to smooth the function
function Xsnew = Smooth(Xs, N, isBoundary, num_iter)
    
    w = 0.3; % relaxation factor
    Xsnew = Xs; % initialize
    %loop through the smoothing function for each iteration
    for n = 1:num_iter
        %loop through every point
        for i = 1:size(Xs,1)
            %ensure we aren't smoothing boundary nodes
            if ~isBoundary(i)
        
                xi = Xs(i,:);
                avg = mean(Xs(N{i},:),1);
        
                Xsnew(i,:) = (1-w)*xi + w*avg;
        
            end
        end
        Xs = Xsnew;
    end
end
