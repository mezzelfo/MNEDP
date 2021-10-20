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
    d = (coords([3,1,2],:)-coords([2,3,1],:)).*[1,-1]; % LUCA
    area = geom.support.TInfo(e).Area;
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
                epsilon/(4*area)*(d(k,2)*d(j,2)+d(k,1)*d(j,1))+...
                1/6*(betax*d(k,2)+betay*d(k,1))+...
                gamma*area/12*(1+k==j);
        end
        center = mean(coords);
        b(jj) = b(jj) + f(center)/3*area;
    end
end
end