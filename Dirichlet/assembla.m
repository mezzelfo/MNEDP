function [A,b] = assembla( ...
    geom, ...
    epsilon, ...
    betax, ...
    betay, ...
    gamma, ...
    f)

ndof = max(geom.pivot.pivot);
A = zeros(ndof,ndof);
b = zeros(ndof,1);

for e=1:geom.nelements.nTriangles
    points_idx = geom.elements.triangles(e,:);
    coords = geom.elements.coordinates(points_idx,:);
    % LUCA: la seguente istruzione calcola dxi,dyi per i=1,2,3 t.c.
    % d(i,1) = dxi e d(i,2) = dyi
    % .*[1,-1] può non servire a patto di ridefinire i dyi come -dyi
    d = (coords([3,1,2],:)-coords([2,3,1],:)).*[1,-1];    
    area = geom.support.TInfo(e).Area;
    center = mean(coords);
    for j=1:3
        jj = geom.pivot.pivot(points_idx(j));
        if jj < 0
            continue
        end
        for k = 1:3
            kk = geom.pivot.pivot(points_idx(k));
            if kk < 0
                continue
            end
            A(jj,kk) = A(jj,kk) + ...
                epsilon(center)/(4*area)*(d(k,2)*d(j,2)+d(k,1)*d(j,1))+...
                1/6*(betax(center)*d(k,2)+betay(center)*d(k,1))+...
                gamma(center)*area/12*(1+k==j);
        end
        b(jj) = b(jj) + f(center)/3*area;
    end
end
end