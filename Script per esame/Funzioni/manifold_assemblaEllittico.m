function [A,b,M] = manifold_assemblaEllittico(Pk,time)

global geom;
global manifold;
global problem;

ndof = max(geom.pivot.pivot);
A = spalloc(ndof,ndof,10*ndof);
if time
    M = spalloc(ndof,ndof,10*ndof);
else
    M = NaN;
end

b = zeros(ndof,1);

% Nodi e pesi per la quadratura nel triangolo di riferimento
[zita,csi,eta,omega] = int_nodes_weights(5);
n_quad = length(omega);

% Funzioni di base e gradienti calcolati nei punti di quadratura
switch Pk
    case 'P1'
        phi = [csi; eta; zita];
        grad_phi = [[1 0];[0 1];[-1 -1]]';
        grad_phi = repmat(grad_phi,1,1,length(omega));
        %     case 'P2'
        %         phi = [2.*csi.*(csi-0.5);2.*eta.*(eta-0.5);2.*zita.*(zita-0.5);4*csi.*eta;4.*eta.*zita;4.*csi.*zita];
        %         grad_phi = zeros(2,6,7);
        %         grad_phi(1,1,:) = -1+4*csi;
        %         grad_phi(2,2,:) = -1+4*eta;
        %         grad_phi(1,3,:) = -3+4*csi+4*eta;
        %         grad_phi(2,3,:) = -3+4*csi+4*eta;
        %         grad_phi(1,4,:) = 4*eta;
        %         grad_phi(2,4,:) = 4*csi;
        %         grad_phi(1,5,:) = -4*eta;
        %         grad_phi(2,5,:) = -4*(-1+csi+2*eta);
        %         grad_phi(1,6,:) = -4*(-1+2*csi+eta);
        %         grad_phi(2,6,:) = -4*csi;
    otherwise
        error('Pk only implemented with P1');
end

Ginv = zeros(2,2,length(csi));
for e = 1:geom.nelements.nTriangles
    points_idx = geom.elements.triangles(e,:);
    coords = geom.elements.coordinates(points_idx,:)';
    d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
    area = geom.support.TInfo(e).Area;
    B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
    invB = inv(B);
    pts = coords(:,end)+B*[csi; eta];
    sqrtdetg = manifold.numeric.sqrtdetg(pts);
    
    for c = 1:n_quad
        Ginv(:,:,c) = manifold.numeric.ginv(pts(:,c));
    end
    epsilon = problem.epsilon(pts);
    beta = problem.beta(pts);
    sigma = problem.sigma(pts);
    f = problem.f(pts);
    BinvTgrad_phi = pagemtimes(invB,'transpose',grad_phi,'none');
    Ginvgrad = pagemtimes(Ginv,BinvTgrad_phi);
    gradGinvgrad = pagemtimes(BinvTgrad_phi,'transpose',Ginvgrad,'none');
    
    betaTgrad = squeeze(sum(reshape(beta,[2,1,7]).*BinvTgrad_phi));
    phiphi = pagemtimes(reshape(phi,[],1,7),reshape(phi,1,[],7));
    
    diffusione = 2*area*squeeze(pagemtimes(omega.*sqrtdetg.*epsilon, permute(gradGinvgrad,[3,1,2])));
    convezione = 2*area*(omega.*sqrtdetg.*betaTgrad) * phi';
    reazione = 2*area*squeeze(pagemtimes(omega.*sqrtdetg.*sigma, permute(phiphi,[3,1,2])));
    total = diffusione+convezione+reazione;
    
    if time
        rho = problem.rho(pts);
        massa = 2*area*squeeze(pagemtimes(omega.*sqrtdetg.*rho, permute(phiphi,[3,1,2])));
    end
    
    for j=1:size(phi,1)
        jj = geom.pivot.pivot(points_idx(j));
        if jj > 0
            for k = 1:size(phi,1)
                kk = geom.pivot.pivot(points_idx(k));
                if kk > 0
                    A(jj,kk) = A(jj,kk) + total(k,j);
                    if time
                        M(jj,kk) = M(jj,kk) + massa(k,j);
                    end
                else % Dirichlet
                    error('No Dirichlet!')
                end
            end
            b(jj) = b(jj) + 2*area*(omega.*sqrtdetg.*f)*phi(j,:)'; % Termine noto
        else % Dirichlet
            error('No Dirichlet!')
        end
    end
end
end