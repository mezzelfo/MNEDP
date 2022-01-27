function u = manifold_assemblaIperbolico(Pk,T,n_passi,beta_Newmark)

global geom;
global manifold;
global problem;

% Fisso il tempo a 0 per calcolare i valori iniziali di f
f_time = problem.f;
problem.f = @(x) f_time(x,0);

[A,f_old_old,M] = manifold_assemblaEllittico(Pk,true);

problem.f = f_time;

delta_t = T/n_passi;

f_old = zeros(size(f_old_old));
f_new = zeros(size(f_old_old));

u = zeros(length(f_new),n_passi+1);

[zita,csi,eta,omega] = int_nodes_weights(5);
phi = [csi; eta; zita];

% Interpolazione della condizione iniziale
for n = 1:geom.nelements.nVertexes
    nn = geom.pivot.pivot(n);
    if nn > 0
        u(nn,1) = problem.iniziale(geom.elements.coordinates(n,:)');
    end
end
% Interpolazione della velocità iniziale
v = zeros(length(f_new),1);
for n = 1:geom.nelements.nVertexes
    nn = geom.pivot.pivot(n);
    if nn > 0
        v(nn) = problem.v_iniziale(geom.elements.coordinates(n,:)');
    end
end

% Genero il primo passo
u(:,2) = u(:,1) + delta_t*v;

% Assemblo f al tempo 1
for e = 1:geom.nelements.nTriangles
    indices = geom.elements.triangles(e,:); % Vertices indices
    coords = geom.elements.coordinates(indices,:)'; % Vertices coordinates
    d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
    area = geom.support.TInfo(e).Area;
    B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
    pts = coords(:,end)+B*[csi; eta];
    
    f = problem.f(pts,delta_t);
    sqrtdetg = manifold.numeric.sqrtdetg(pts);
    ii = geom.pivot.pivot(indices);
    f_old(ii(ii > 0)) = f_old(ii(ii > 0)) + 2*area*sum(omega.*sqrtdetg.*f.*phi(ii > 0,:),2);
end

m_Newmark = M + beta_Newmark*delta_t^2*A;
m_Newmark2 = 2*M - (1-2*beta_Newmark)*delta_t^2*A;

% Avanzamento in tempo con il metodo di Newmark
for t = 1:n_passi-1
    %     for e = 1:geom.nelements.nTriangles
    %         indices = geom.elements.triangles(e,:); % Vertices indices
    %         coords = geom.elements.coordinates(indices,:)'; % Vertices coordinates
    %         d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
    %         area = geom.support.TInfo(e).Area;
    %         B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
    %         pts = coords(:,end)+B*[csi; eta];
    %
    %         f = problem.f(pts,t*delta_t);
    %         sqrtdetg = manifold.numeric.sqrtdetg(pts);
    %         ii = geom.pivot.pivot(indices);
    %         f_new(ii(ii > 0)) = f_new(ii(ii > 0)) + 2*area*sum(omega.*sqrtdetg.*f.*phi(ii > 0,:),2);
    %     end
    %     u(:,t+2) = m_Newmark\(m_Newmark2*u(:,t+1) - m_Newmark*u(:,t) + ...
    %         delta_t^2*(beta_Newmark*f_new + (1-2*beta_Newmark)*f_old + beta_Newmark*f_old_old));
    %     f_old_old = f_old;
    %     f_old = f_new;
    %     f_new = 0*f_new;
    
    %SPEED: se f è zero
    u(:,t+2) = m_Newmark\(m_Newmark2*u(:,t+1) - m_Newmark*u(:,t));
end
end