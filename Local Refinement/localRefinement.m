function [nodes,NE] = localRefinement(nodes,NE,BE,IE,X_BE_spline,Title)

    smoothIter = 5;
    
    done = false;
    
    while ~done
    
        nElem = size(NE,1);
        nNode = size(nodes,1);
    
        %% ============================================================
        % Flag edges using sizing function
        %% ============================================================
    
        [Iflagged,Bflagged] = FlagEdge(nodes,IE,BE,X_BE_spline);
    
        if ~any(Iflagged) && ~any(Bflagged)
            break;
        end
    
        %% ============================================================
        % Element flags
        %% ============================================================
    
        elemFlag = false(nElem,1);
    
        for i=find(Iflagged)'
            elemFlag(IE(i,3)) = true;
            elemFlag(IE(i,4)) = true;
        end
    
        for i=find(Bflagged)'
            elemFlag(BE(i,3)) = true;
        end
    
        %% ============================================================
        % Edge midpoint table
        %% ============================================================
    
        edgeMid = containers.Map('KeyType','char','ValueType','double');
    
        newNodes = nodes;
        newElem = [];
    
        %% ============================================================
        % Refinement loop
        %% ============================================================
    
        for e=1:nElem
    
            tri = NE(e,:);
    
            if ~elemFlag(e)
                newElem(end+1,:) = tri;
                continue;
            end
    
            % edges of triangle
            E = [tri([1 2]);
                 tri([2 3]);
                 tri([3 1])];
    
            mids = zeros(1,3);
            split = false(1,3);
    
            for k=1:3
                n1=E(k,1); n2=E(k,2);
                key = sprintf('%d_%d',min(n1,n2),max(n1,n2));
    
                hcur = norm(nodes(n2,:)-nodes(n1,:));
                hdes = sizing_function.sizing((nodes(n1,:)+nodes(n2,:))/2,X_BE_spline);
    
                if hcur > hdes
                    split(k)=true;
    
                    if ~isKey(edgeMid,key)
                        newNodes(end+1,:)=(nodes(n1,:)+nodes(n2,:))/2;
                        edgeMid(key)=size(newNodes,1);
                    end
    
                    mids(k)=edgeMid(key);
                end
            end
    
            ns=sum(split);
    
            %% =====================================================
            % 3-edge split
            %% =====================================================
    
            if ns==3
                a=mids(1); b=mids(2); c=mids(3);
    
                newElem=[newElem;
                    tri(1) a c;
                    a tri(2) b;
                    c b tri(3);
                    a b c];
    
            %% =====================================================
            % 2-edge split
            %% =====================================================
    
            elseif ns==2
    
                idx=find(split);
                m1=mids(idx(1));
                m2=mids(idx(2));
    
                opp=setdiff(1:3,idx);
                v=tri(opp);
    
                newElem=[newElem;
                    v m1 m2;
                    E(idx(1),1) m1 v;
                    E(idx(2),2) m2 v];
    
            %% =====================================================
            % 1-edge split (longest angle rule)
            %% =====================================================
    
            elseif ns==1
    
                k=find(split);
                m=mids(k);
    
                verts=[tri(k) tri(mod(k,3)+1)];
                opp=tri(mod(k+1,3)+1);
    
                newElem=[newElem;
                    verts(1) m opp;
                    m verts(2) opp];
    
            else
                newElem(end+1,:)=tri;
            end
        end
    
        %% ============================================================
        % Update mesh
        %% ============================================================
    
        nodes=newNodes;
        NE=newElem;
    
        %% ============================================================
        % Smoothing
        %% ============================================================
    
        N=Connectivity(NE,size(nodes,1));
        isBoundary=BuildBoundaryNodes(BE,NE,size(nodes,1));
        nodes=Smooth(nodes,N,isBoundary,smoothIter);
    
        %% ============================================================
        % Rebuild edges for next iteration (user routine)
        %% ============================================================
    
        [IE,BE]=BuildEdges(NE);
    
    end

end