% performs uniform refinement N number of times
function uniformRefinement(grifile,nRef)

    newgrifile = [grifile(1:end-4),num2str(nRef),grifile(end-3:end)];

    % loop refinement n times
    for ref = 1:nRef
        %% load gri  grifile = 'test.gri';
        fid = fopen(grifile, 'r');
        A = fscanf(fid,'%d', 3);
        nNode    = A(1);
        nElemTot = A(2);
        Dim      = A(3);

        % read through nodes
        nodes = zeros(nNode, Dim);
        for iNode = 1:nNode
            A = fscanf(fid, '%lf', 2);
            nodes(iNode,:) = A(1:2)';
        end

        % read through boundary info
        A = fscanf(fid, '%d', 1);
        nBGroup = A(1);
        edgeTypes = sparse(nNode,nNode);
        for iBGroup = 1:nBGroup
            fgets(fid);
            sLine = fgets(fid);
            [nBFace(iBGroup), nf(iBGroup), Title(iBGroup)] = strread(sLine, '%d %d %s');
            for iBFace = 1:nBFace(iBGroup)
                A = fscanf(fid, '%d', nf(iBGroup));
                edgeTypes(A(1),A(2)) = iBGroup;
                edgeTypes(A(2),A(1)) = iBGroup;
            end
        end

        % read through elements
        curtot = 0;
        NE = zeros(nElemTot,3);
        nElemGroup = 0;
        while (curtot ~= nElemTot)
            nElemGroup = nElemGroup+1;
            fgets(fid);
            sLine = fgets(fid);
            [nElem(nElemGroup), Order(nElemGroup), Basis(nElemGroup)] = strread(sLine, '%d %d %s');
            % save these
            switch Basis{1}
                case 'TriLagrange'
                    nnode = (Order+1)*(Order+2)/2;
                    nedge = 3;
                otherwise
                    error('element type not understood');
            end
            for elem = 1:nElem
                A = fscanf(fid, '%d', nnode);
                for edge=1:nedge
                    NE(elem,:) = A;
                end
            end
            curtot = curtot + nElem;
        end

        % read through periodic groups
        fgets(fid);
        sLine = fgets(fid);
        [nPG, PeriodicGroup] = strread(sLine, '%d %s');
        for iPGroups = 1:nPG
            sLine = fgets(fid);
            [nPGNode(iPGroups), Periodicity(iPGroups)] = strread(sLine, '%d %s');
            for iPFaces = 1:nPGNode
                A = fscanf(fid, '%d', 2);
                NP(iPFaces,:,iPGroups) = A;
            end
            fgets(fid); % needs to bere here for multiple periodic groups
        end


        %% build new mesh
        visitedEdges = sparse(nNode,nNode);
        edgerec = zeros(1,2);
        nNodeNew = nNode;
        newNodes = nodes;
        newNodeNums = zeros(nedge,1);
        nElemNew = nElemTot*4;
        newElem = zeros(nElemNew,3);
        newEdges = zeros(max(nBFace)*2,2,nBGroup);
        edgeIndex = ones(nBGroup,1);
        periodicNodes = sparse(nNode,nNode);
        i = 1;
        for elem = 1:nElemTot
            elemnodes = [NE(elem,:),NE(elem,1)];
            % on each edge see if edge has been visited before with sparce,
            %   if not, incriment total node number and save to sparce
            for edge = 1:nedge
                edgerec(1) = elemnodes(edge);
                edgerec(2) = elemnodes(edge+1);
                if visitedEdges(edgerec(1),edgerec(2))==0
                    nNodeNew = nNodeNew + 1;
                    visitedEdges(edgerec(1),edgerec(2)) = nNodeNew;
                    visitedEdges(edgerec(2),edgerec(1)) = nNodeNew;
                    % edges
                    tmpBGroup = edgeTypes(edgerec(1),edgerec(2));
                    if tmpBGroup ~= 0
                        newEdges(edgeIndex(tmpBGroup),:,tmpBGroup) = [edgerec(1),nNodeNew];
                        newEdges(edgeIndex(tmpBGroup)+1,:,tmpBGroup) = [nNodeNew,edgerec(2)];
                        edgeIndex(tmpBGroup) = edgeIndex(tmpBGroup)+2;
                        % periodic stuff
                        for iPGroups = 1:nPG
                            inPGroup = (sum(NP(:,:,iPGroups)==edgerec(1)))+(sum(NP(:,:,iPGroups)==edgerec(2)));
                            if sum(inPGroup==[2,2])
                                [row(1),col] = find(NP(:,:,iPGroups)==edgerec(1));
                                [row(2),col] = find(NP(:,:,iPGroups)==edgerec(2));
                                pairedNodes = [NP(row(1),rem(col,2)+1,iPGroups),NP(row(2),rem(col,2)+1,iPGroups)];
                                if periodicNodes(pairedNodes(1),pairedNodes(2)) == 0
                                    nPGNode(iPGroups) = nPGNode(iPGroups)+1;
                                    periodicNodes(edgerec(1),edgerec(2)) = nPGNode(iPGroups);
                                    periodicNodes(edgerec(2),edgerec(1)) = nPGNode(iPGroups);
                                    NP(nPGNode(iPGroups),col,iPGroups) = nNodeNew;
                                else
                                    NP(periodicNodes(pairedNodes(1),pairedNodes(2)),col,iPGroups) = nNodeNew;
                                end
                            end
                        end
                        % THIS CALCS NEW XY if an edge, so add projection here
                        newNodes(nNodeNew,:) = (nodes(elemnodes(edge),:)+nodes(elemnodes(edge+1),:))/2;
                    else
                        %if not an edge, jsut avg x & y values
                        newNodes(nNodeNew,:) = (nodes(elemnodes(edge),:)+nodes(elemnodes(edge+1),:))/2;
                    end
                    newNodeNums(edge) = nNodeNew;
                else
                    newNodeNums(edge) = visitedEdges(edgerec(1),edgerec(2));
                end
            end
            % ensures lowest node# will be first
            [~,I] = min(newNodeNums);
            temp = [1,2,3,1,2];
            % adds to new element list
            newElem((elem*4-3),:,i) = [elemnodes(1),newNodeNums(1),newNodeNums(3)];
            newElem((elem*4-2),:,i) = [elemnodes(2),newNodeNums(2),newNodeNums(1)];
            newElem((elem*4-1),:,i) = [newNodeNums(temp(I)),newNodeNums(temp(I+1)),newNodeNums(temp(I+2))];
            newElem((elem*4-0),:,i) = [elemnodes(3),newNodeNums(3),newNodeNums(2)];
        end

        
        %% set old to new
        nNode = nNodeNew;
        nElemTot = nElemNew;
        nodes = newNodes;
        nBFace = nBFace*2;
        NB = newEdges;
        nElem = nElem*4;
        NE = newElem;


        %% Save to new .gri file
        % open file for writing
        fname = sprintf(newgrifile);
        fid = fopen(fname, 'w');

        % top
        fprintf(fid, '%d %d %d\n', nNode, nElemTot, Dim);

        % nodes
        for i = 1:nNode
            fprintf(fid, '%20.15f %20.15f\n', nodes(i,1), nodes(i,2));
        end

        % boundaries
        fprintf(fid, '%d\n', nBGroup);
        for i = 1:nBGroup
            fprintf(fid, '%d %d %s\n', nBFace(i), nf(i), Title{1,i});
            for j = 1:nBFace(i)
                fprintf(fid, '%d %d\n', NB(j,1,i), NB(j,2,i));
            end
        end

        % elements
        for i = 1:nElemGroup
            fprintf(fid, '%d %d %s\n', nElem(i), Order(i), Basis{1,i});
            for j = 1:nElem(i)
                fprintf(fid, '%d %d %d\n', NE(j,1,i), NE(j,2,i), NE(j,3,i));
            end
        end

        % periodics
        fprintf(fid, '%d %s\n', nPG, PeriodicGroup{1,1});
        for i = 1:nPG
            % nPGNode(i) Periodicity(i)
            fprintf(fid, '%d %s\n', nPGNode(i), Periodicity{1,i});
            for j = 1:nPGNode(i)
                % NP(i,j,1) NP(i,j,2)
                fprintf(fid, '%d %d \n', NP(j,1,i), NP(j,2,i));
            end
        end
        
        % set new to old file
        grifile = newgrifile;
    end
end
