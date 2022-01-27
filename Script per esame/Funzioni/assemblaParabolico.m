function [u,uD] = assemblaParabolico(Pk,T,n_passi,SUPG,MassLumping,SpeedMode)

global geom;
global problem;

% Fisso il tempo a 0 per calcolare i valori iniziali di f,uD e fN
f_time = problem.f;
bordo_dirichlet_time = problem.bordo_dirichlet;
bordo_neumann_time = problem.bordo_neumann;

problem.f = @(x) f_time(x,0);
problem.bordo_dirichlet = @(x,marker) bordo_dirichlet_time(x,0,marker);
problem.bordo_neumann = @(x,marker) bordo_neumann_time(x,0,marker);

[A,f_old,AD,~,fN_old,M,MD] = assemblaEllittico(Pk,SUPG,MassLumping,true);
problem.f = f_time;
problem.bordo_dirichlet = bordo_dirichlet_time;
problem.bordo_neumann = bordo_neumann_time;

delta_t = T/n_passi;

f_new = zeros(size(f_old));
fN_new = zeros(size(fN_old));

u = zeros(length(f_new),n_passi+1); % Soluzione sui pivot nel tempo
uD = zeros(size(AD,2),n_passi+1); % Soluzione al bordo nel tempo

% Pesi e nodi per formula quadratura 2D
[zita,csi,eta,omega] = int_nodes_weights(5);

% Pesi e nodi per formula quadratura 1D
csi_1D = [-sqrt(5+2*sqrt(10/7))/3, -sqrt(5-2*sqrt(10/7))/3, 0, sqrt(5-2*sqrt(10/7))/3, sqrt(5+2*sqrt(10/7))/3];
omega_1D = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900];
csi_1D = csi_1D/2+0.5;
omega_1D = omega_1D/2;

switch Pk
    case 'P1'
        phi = [csi; eta; zita];
        P_1D = [1-csi_1D; csi_1D];
    case 'P2'
        phi = [2.*csi.*(csi-0.5);2.*eta.*(eta-0.5);2.*zita.*(zita-0.5);4*csi.*eta;4.*eta.*zita;4.*csi.*zita];
        P_1D = [2.*(csi_1D-1).*(csi_1D-0.5); 2.*csi_1D.*(csi_1D-0.5); -4.*csi_1D.*(csi_1D-1)];
    otherwise
        error('Funzione implementata solo per P1 e P2');
end

% Interpolazione della condizione iniziale
u(:,1) = problem.iniziale(geom.elements.coordinates(geom.pivot.pivot>0,:)');
uD(:,1) = problem.iniziale(geom.elements.coordinates(geom.pivot.pivot<0,:)');

% Calcolo delle matrici utilizzate nell'avanzamento in tempo
Kp = M + A*delta_t/2;
Km = M - A*delta_t/2;
KpD = MD + AD*delta_t/2;
KmD = MD - AD*delta_t/2;
% Indici per il calcolo delle condizioni di Neumann
dof_border_position = logical(1:size(geom.elements.borders,2));
dof_border_position([3,4]) = 0;


perm = symrcm(Kp);

Kp = Kp(perm,perm);
Km = Km(perm,perm);
KpD = KpD(perm,:);
KmD = KmD(perm,:);

% Avanzamento in tempo con il metodo di Crank-Nicolson
for t = 1:n_passi
    if SpeedMode == true
        u(perm,t+1) = Kp\(Km*u(perm,t));
        %u(:,t+1) = Kp\(Km*u(:,t));
    else
        for e = 1:geom.nelements.nTriangles % Aggiorna la forzante
            indices = geom.elements.triangles(e,:); % Vertices indices
            coords = geom.elements.coordinates(indices,:)'; % Vertices coordinates
            d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
            area = geom.support.TInfo(e).Area;
            B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
            pts = coords(:,3)+B*[csi; eta];
            
            ii = geom.pivot.pivot(indices);
            
            fpts = problem.f(pts,t*delta_t);
            
            f_new(ii(ii > 0)) = f_new(ii(ii > 0)) + (2*area*omega*(fpts.*phi(ii > 0,:))')';
        end
        for j = 1:length(geom.pivot.Di) % Aggiorna le condizioni di dirichlet
            j = geom.pivot.Di(j,:);
            jj = -geom.pivot.pivot(j(1));
            coords = geom.elements.coordinates(j(1),:)';
            uD(jj,t+1) =  problem.bordo_dirichlet(coords,t*delta_t, j(2));
        end
        for e = 1:length(geom.pivot.Ne) % Aggiorna le condizioni di Neumann
            e = geom.pivot.Ne(e,:);
            points_idx = geom.elements.borders(e(1),dof_border_position); %[start node, end node, middle node]
            coords = geom.elements.coordinates(points_idx,:)';
            
            border_length = norm(coords(:,1)-coords(:,2));
            
            %from csi_1D to actual points along border and then get neumann values
            pts = coords(:,1) + (coords(:,2)-coords(:,1))*csi_1D;
            gamma_Ne = problem.bordo_neumann(pts,t*delta_t, e(2)); %e(2) is marker
            
            % Numerical integration of the three basis functions
            contribs = border_length * (omega_1D .* gamma_Ne) * P_1D';
            
            jj = geom.pivot.pivot(points_idx);
            fN_new(jj(jj > 0)) = fN_new(jj(jj > 0)) + contribs(jj > 0)';
        end
        %u(:,t+1) = Kp\(Km*u(:,t) - KpD*uD(:,t+1) + KmD*uD(:,t) + delta_t/2*(f_old + f_new + fN_old + fN_new));
        f_tmp = f_old + f_new + fN_old + fN_new;
        u(perm,t+1) = Kp\(Km*u(perm,t) - KpD*uD(:,t+1) + KmD*uD(:,t) + delta_t/2*f_tmp(perm));
        f_old = f_new;
        f_new = zeros(size(f_new));
        fN_old = fN_new;
        fN_new = zeros(size(fN_new));
    end
end
end