function [A,b,A_dirichlet,u_dirichlet,b_neumann] = assembla( ...
    epsilon, ...
    betax, ...
    betay, ...
    gamma, ...
    f,...
    bordo_dirichlet,...
    bordo_neumann)

global geom

ndof = max(geom.pivot.pivot);

A = spalloc(ndof,ndof, 10*ndof);
b = zeros(ndof,1);

ndirichlet = -min(geom.pivot.pivot);
A_dirichlet = spalloc(ndof,ndirichlet, 10*ndof);
u_dirichlet = zeros(ndirichlet,1);

b_neumann = zeros(ndof,1);

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
        if jj > 0
            for k = 1:3
                kk = geom.pivot.pivot(points_idx(k));
                if kk > 0
                    A(jj,kk) = A(jj,kk) + ...
                        epsilon(center)/(4*area)*(d(k,2)*d(j,2)+d(k,1)*d(j,1))+...
                        1/6*(betax(center)*d(k,2)+betay(center)*d(k,1))+...
                        gamma(center)*area/12*(1+k==j);
                else
                    A_dirichlet(jj,-kk) = A_dirichlet(jj,-kk) + ...
                        epsilon(center)/(4*area)*(d(k,2)*d(j,2)+d(k,1)*d(j,1))+...
                        1/6*(betax(center)*d(k,2)+betay(center)*d(k,1))+...
                        gamma(center)*area/12*(1+k==j);
                end
            end
            b(jj) = b(jj) + f(center)/3*area;
        end
    end
end

for j = 1:length(geom.pivot.Di)
    j = geom.pivot.Di(j,:);
    jj = -geom.pivot.pivot(j(1));
    coords = geom.elements.coordinates(j(1));
    marker = j(2);
    u_dirichlet(jj) = bordo_dirichlet(coords,marker);
end

for e = 1:length(geom.pivot.Ne)
    e = geom.pivot.Ne(e,:);
    points_idx = geom.elements.borders(e(1),[1,2]);
    coords = geom.elements.coordinates(points_idx,:);
    lunghezza = sqrt(sum((coords(1,:)-coords(2,:)).^2));
    marker = e(2);
    gn = bordo_neumann(coords, marker);
    contribs = lunghezza * [[2, 1];[1, 2]]/6 * gn;
    jj = geom.pivot.pivot(points_idx);
    b_neumann(jj(jj > 0)) = b_neumann(jj(jj > 0)) + contribs(jj > 0);   
end


end