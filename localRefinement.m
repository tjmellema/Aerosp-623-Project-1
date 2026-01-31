function [nodes,NE] = localRefinement(nodes,NE,BE,IE,...
                       X_BE_spline_top,X_BE_spline_bottom,...
                       periodic_pairs)

smoothIter = 5;

% ===============================
% Build periodic node map
% ===============================

pmap = containers.Map('KeyType','double','ValueType','double');

for i=1:size(periodic_pairs,1)
    a = periodic_pairs(i,1);
    b = periodic_pairs(i,2);
    pmap(a)=b;
    pmap(b)=a;
end

% ===============================
% Track boundary nodes explicitly
% ===============================

isBoundaryNode = false(size(nodes,1),1);
isBoundaryNode(unique(BE(:,1:2))) = true;

done = false;

while ~done

    nElem = size(NE,1);

    % =====================================
    % Flag edges
    % =====================================

    [Iflagged,Bflagged] = FlagEdge(nodes,IE,BE,...
                          X_BE_spline_top,X_BE_spline_bottom);

    % =====================================
    % Enforce periodic matching on BE
    % =====================================

    for i=1:size(BE,1)

        n1=BE(i,1); n2=BE(i,2);

        if isKey(pmap,n1) && isKey(pmap,n2)

            pn1=pmap(n1); pn2=pmap(n2);

            for j=1:size(BE,1)
                if (BE(j,1)==pn1 && BE(j,2)==pn2) || ...
                   (BE(j,1)==pn2 && BE(j,2)==pn1)

                    if Bflagged(i), Bflagged(j)=true; end
                    if Bflagged(j), Bflagged(i)=true; end
                end
            end
        end
    end

    if ~any(Iflagged) && ~any(Bflagged)
        break;
    end

    % =====================================
    % Element flags
    % =====================================

    elemFlag=false(nElem,1);

    for i=find(Iflagged)'
        elemFlag(IE(i,3))=true;
        elemFlag(IE(i,4))=true;
    end

    for i=find(Bflagged)'
        elemFlag(BE(i,3))=true;
    end

    % =====================================
    % Midpoint table
    % =====================================

    edgeMid=containers.Map('KeyType','char','ValueType','double');

    newNodes=nodes;
    newElem=[];
    newBoundary = isBoundaryNode;

    % =====================================
    % Refinement loop
    % =====================================

    for e=1:nElem

        tri=NE(e,:);

        if ~elemFlag(e)
            newElem(end+1,:)=tri;
            continue
        end

        E=[tri([1 2]);
           tri([2 3]);
           tri([3 1])];

        mids=zeros(1,3);
        split=false(1,3);

        for k=1:3

            n1=E(k,1); n2=E(k,2);

            key=sprintf('%d_%d',min(n1,n2),max(n1,n2));

            hcur=norm(nodes(n2,:)-nodes(n1,:));
            hdes=sizing(0.5*(nodes(n1,:)+nodes(n2,:)),...
                         X_BE_spline_top,X_BE_spline_bottom);

            if hcur>hdes

                split(k)=true;

                if ~isKey(edgeMid,key)

                    mid = 0.5*(nodes(n1,:)+nodes(n2,:));
                    newNodes(end+1,:) = mid;
                    edgeMid(key) = size(newNodes,1);

                    % --- boundary tagging ---
                    isB = isBoundaryNode(n1) && isBoundaryNode(n2);
                    newBoundary(end+1) = isB;

                    % --- periodic midpoint creation ---
                    if isKey(pmap,n1) && isKey(pmap,n2)

                        pn1=pmap(n1); pn2=pmap(n2);
                        pkey=sprintf('%d_%d',min(pn1,pn2),max(pn1,pn2));

                        pmid = 0.5*(nodes(pn1,:)+nodes(pn2,:));
                        newNodes(end+1,:) = pmid;
                        edgeMid(pkey) = size(newNodes,1);

                        % periodic midpoints always boundary
                        newBoundary(end+1) = true;
                    end
                end

                mids(k)=edgeMid(key);
            end
        end

        ns=sum(split);

        if ns==3
            a=mids(1); b=mids(2); c=mids(3);
            newElem=[newElem;
                tri(1) a c;
                a tri(2) b;
                c b tri(3);
                a b c];

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

    nodes = newNodes;
    NE    = newElem;
    isBoundaryNode = newBoundary;

    % =====================================
    % Smoothing (boundary excluded)
    % =====================================

    N = Connectivity(NE,size(nodes,1));
    nodes = Smooth(nodes,N,isBoundaryNode,smoothIter);

    % rebuild edges
    [IE,BE] = BuildEdges(NE);

end

end
