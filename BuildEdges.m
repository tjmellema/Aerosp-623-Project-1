function [IE,BE] = BuildEdges(NE)
% BUILD_EDGES
%
% From element connectivity NE (nElem x 3),
% construct:
%
% IE(i,:) = [n1 n2 elemL elemR]
% BE(i,:) = [n1 n2 elem]
%
% Author: ChatGPT

nElem = size(NE,1);

% map: "n1_n2" → [elem, localEdge]
edgeMap = containers.Map('KeyType','char','ValueType','any');

IE = [];
BE = [];

% local edges of triangle
loc = [1 2;
       2 3;
       3 1];

for e = 1:nElem

    tri = NE(e,:);

    for k = 1:3

        n1 = tri(loc(k,1));
        n2 = tri(loc(k,2));

        % enforce ordering
        a = min(n1,n2);
        b = max(n1,n2);

        key = sprintf('%d_%d',a,b);

        if ~isKey(edgeMap,key)
            edgeMap(key) = [e];
        else
            edgeMap(key) = [edgeMap(key), e];
        end

    end
end

keysE = keys(edgeMap);

for i=1:length(keysE)

    key = keysE{i};
    elems = edgeMap(key);

    tok = sscanf(key,'%d_%d');

    if length(elems)==2
        IE(end+1,:) = [tok(1) tok(2) elems(1) elems(2)];
    else
        BE(end+1,:) = [tok(1) tok(2) elems(1)];
    end
end

end