function prepP2()
global geom

nVertexes = geom.nelements.nVertexes;

for b = 1:geom.nelements.nBorders
    nVertexes = nVertexes+1;
    idx_vertices = geom.elements.borders(b,1:2);
    idx_elements = geom.elements.borders(b,3:4);

    %Add new node to adjacent triangles
    for element = idx_elements
        if(element ~= -1)
            pos = geom.elements.triangles(element,1:3) == idx_vertices';
            pos = sum(pos * [5;3;7]/2);
            assert((pos == 4) | (pos == 5) | (pos == 6));
            geom.elements.triangles(element,pos) = nVertexes;
        end
    end

    %Get correct marker for the new node
    edge_marker = geom.support.BInfo(b,3);
    if mod(edge_marker, 2) == 0
        %Neumann condition or no border condition at all
        geom.pivot.nodelist(nVertexes) = 0;
        geom.pivot.pivot(nVertexes) = max(geom.pivot.pivot)+1;
    else
        %Dirichlet condition
        geom.pivot.nodelist(nVertexes)= edge_marker;
        geom.pivot.Di(end+1,:) = [nVertexes, edge_marker];
        geom.pivot.pivot(nVertexes) = min(geom.pivot.pivot)-1;
    end
    
    %TODO: add new coord in geom.elements.coordinates
end
end