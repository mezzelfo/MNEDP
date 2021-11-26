function [A,b,A_dirichlet,u_dirichlet,b_neumann] = assembla( ...
    epsilon, ...
    beta, ...
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

[zita,csi,eta,omega] = int_nodes_weights(5); %TODO: da rendere globale
% Calcoliamo una volta per tutte le funzioni di base su tutti i punti
% della quadratura sul triangolo di riferimento

%%P1
%P = [csi; eta; zita]';
%grad_P = [[1 0];[0 1];[-1 -1]]';
%grad_P = repmat(grad_P,1,1,length(omega));

%%P2
P = [2.*csi.*(csi-0.5);2.*eta.*(eta-0.5);2.*zita.*(zita-0.5);4*csi.*eta;4.*eta.*zita;4.*csi.*zita]';
grad_P = zeros(2,6,7);
grad_P(1,1,:) = -1+4*csi;
grad_P(2,2,:) = -1+4*eta;
grad_P(1,3,:) = -3+4*csi+4*eta;
grad_P(2,3,:) = -3+4*csi+4*eta;
grad_P(1,4,:) = 4*eta;
grad_P(2,4,:) = 4*csi;
grad_P(1,5,:) = -4*eta;
grad_P(2,5,:) = -4*(-1+csi+2*eta);
grad_P(1,6,:) = -4*(-1+2*csi+eta);
grad_P(2,6,:) = -4*csi;

csi_1D = [-sqrt(5+2*sqrt(10/7))/3, -sqrt(5-2*sqrt(10/7))/3, 0, sqrt(5-2*sqrt(10/7))/3, sqrt(5+2*sqrt(10/7))/3];
omega_1D = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900];

csi_1D = csi_1D/2+0.5;
omega_1D = omega_1D/2;

P_1D = [2.*(csi_1D-1).*(csi_1D-0.5); 2.*csi_1D.*(csi_1D-0.5); -4.*csi_1D.*(csi_1D-1)]';

for e=1:geom.nelements.nTriangles
    dof = geom.elements.triangles(e,:);
    points_idx = geom.elements.triangles(e,1:3);
    coords = geom.elements.coordinates(points_idx,:)';
    d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
    area = geom.support.TInfo(e).Area;
    center = mean(coords,2);
    
    B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
    BinvTgrad_P = pagemtimes(inv(B)',grad_P);
    gradgrad = pagemtimes(pagetranspose(BinvTgrad_P),BinvTgrad_P);
    
    pts = coords(:,end)+B*[csi; eta];
    
    betaTgrad = squeeze(sum(reshape(beta(pts),[2,1,7]).*BinvTgrad_P));
    
    mk = 1/3;
    h = sqrt(max(sum(d.^2,1)));
    Pe = mk*norm(beta(center))*h/(2*epsilon(center));
    tau = h/(2*norm(beta(center)));
    if Pe < 1
        tau = mk*h^2/(4*epsilon(center));
    end
    
    tau = 0; %TODO: SUPG sembra sbagliato
    
    for j=1:size(P,2)
        jj = geom.pivot.pivot(dof(j));
        if jj > 0
            for k = 1:size(P,2)
                contrib_diff = 2*area*(omega .* epsilon(pts)) * squeeze(gradgrad(k,j,:));
                contrib_conv = 2*area*(omega .* squeeze(betaTgrad(k,:))) * P(:,j);
                contrib_SUPG = 2*area*tau*omega * (betaTgrad(k,:) .* betaTgrad(j,:))';
                contrib_SUPG = contrib_SUPG + 0; %TODO aggiungere ordine 2
                kk = geom.pivot.pivot(dof(k));
                if kk > 0
                    A(jj,kk) = A(jj,kk) + contrib_diff + contrib_conv + contrib_SUPG;
                else
                    A_dirichlet(jj,-kk) = A_dirichlet(jj,-kk) + contrib_diff + contrib_conv + contrib_SUPG;
                end
            end
            b(jj) = b(jj) + ...
                2*area*(omega.*f(pts))*P(:,j) + ...
                2*area*tau*(omega.*f(pts))*betaTgrad(j,:)';
        end
    end
end

for j = 1:length(geom.pivot.Di)
    j = geom.pivot.Di(j,:);
    jj = -geom.pivot.pivot(j(1));
    coords = geom.elements.coordinates(j(1));
    u_dirichlet(jj) = bordo_dirichlet(coords,j(2)); %j(2) is marker
end


for e = 1:length(geom.pivot.Ne)
    e = geom.pivot.Ne(e,:);
    points_idx = geom.elements.borders(e(1),[1,2,5]); %[start node, end node, middle node]
    coords = geom.elements.coordinates(points_idx,:)';
    
    border_length = norm(coords(:,1)-coords(:,2));
    
    %from csi_1D to actual points along border and then get neumann values
    pts = coords(:,1) + (coords(:,2)-coords(:,1))*csi_1D;
    gamma_Ne = bordo_neumann(pts, e(2)); %e(2) is marker
    
    % Numerical integration of the three basis functions
    contribs = border_length * (omega_1D .* gamma_Ne) * P_1D;
    
    % Accumulate
    jj = geom.pivot.pivot(points_idx);
    b_neumann(jj(jj > 0)) = b_neumann(jj(jj > 0)) + contribs(jj > 0)';
end


end