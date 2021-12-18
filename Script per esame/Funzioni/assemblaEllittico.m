function [A,b,A_dirichlet,u_dirichlet,b_neumann,mat_massa,mat_massa_dirichlet] = assemblaEllittico(Pk,SUPG,MassLumping,massa)
global problem
global geom

ndof = max(geom.pivot.pivot);
ndirichlet = -min(geom.pivot.pivot);

A = spalloc(ndof,ndof, 10*ndof);
b = zeros(ndof,1);
mat_massa = spalloc(ndof,ndof,10*ndof);

if ~massa
    problem.rho = @(x) 0;
end

A_dirichlet = spalloc(ndof,ndirichlet, 10*ndof);
u_dirichlet = zeros(ndirichlet,1);
mat_massa_dirichlet = spalloc(ndof,ndirichlet,10*ndof);

b_neumann = zeros(ndof,1);

%pesi e nodi per formula quadratura 2D
[zita,csi,eta,omega] = int_nodes_weights(5);

%pesi e nodi per formula quadratura 1D
csi_1D = [-sqrt(5+2*sqrt(10/7))/3, -sqrt(5-2*sqrt(10/7))/3, 0, sqrt(5-2*sqrt(10/7))/3, sqrt(5+2*sqrt(10/7))/3];
omega_1D = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900];
csi_1D = csi_1D/2+0.5;
omega_1D = omega_1D/2;

switch Pk
    case 'P1'
        phi = [csi; eta; zita];
        grad_phi = [[1 0];[0 1];[-1 -1]]';
        grad_phi = repmat(grad_phi,1,1,length(omega));
        hessian_phi = zeros(2,2,3);
        P_1D = [1-csi_1D; csi_1D];
        mk = 1/3; % per SUPG
    case 'P2'
        phi = [2.*csi.*(csi-0.5);2.*eta.*(eta-0.5);2.*zita.*(zita-0.5);4*csi.*eta;4.*eta.*zita;4.*csi.*zita];
        grad_phi = zeros(2,6,7);
        grad_phi(1,1,:) = -1+4*csi;
        grad_phi(2,2,:) = -1+4*eta;
        grad_phi(1,3,:) = -3+4*csi+4*eta;
        grad_phi(2,3,:) = -3+4*csi+4*eta;
        grad_phi(1,4,:) = 4*eta;
        grad_phi(2,4,:) = 4*csi;
        grad_phi(1,5,:) = -4*eta;
        grad_phi(2,5,:) = -4*(-1+csi+2*eta);
        grad_phi(1,6,:) = -4*(-1+2*csi+eta);
        grad_phi(2,6,:) = -4*csi;
        hessian_phi = zeros(2,2,6); %Constant so no need for fourth dimension
        hessian_phi(:,:,1) = [4 0;0 0];
        hessian_phi(:,:,2) = [0 0;0 4];
        hessian_phi(:,:,3) = [4 4;4 4];
        hessian_phi(:,:,4) = [0 4;4 0];
        hessian_phi(:,:,5) = [0 -4;-4 -8];
        hessian_phi(:,:,6) = [-8 -4;-4 0];
        P_1D = [2.*(csi_1D-1).*(csi_1D-0.5); 2.*csi_1D.*(csi_1D-0.5); -4.*csi_1D.*(csi_1D-1)];
        mk = 1/24; % per SUPG
    otherwise
        error('Funzione implementata solo per P1 e P2');
end


for e=1:geom.nelements.nTriangles
    dof = geom.elements.triangles(e,:);
    points_idx = geom.elements.triangles(e,1:3); %TODO can be eliminated
    coords = geom.elements.coordinates(points_idx,:)';
    d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
    area = geom.support.TInfo(e).Area;
    center = mean(coords,2);
    
    B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
    invB = inv(B);
    
    pts = coords(:,end)+B*[csi; eta];
    epsilonpts = problem.epsilon(pts);
    betapts = problem.beta(pts);
    sigmapts = problem.sigma(pts);
    fpts = problem.f(pts);
    

    BinvTgrad_phi = pagemtimes(invB',grad_phi);
    gradgrad = pagemtimes(pagetranspose(BinvTgrad_phi),BinvTgrad_phi);
    betaTgrad = squeeze(sum(reshape(betapts,[2,1,7]).*BinvTgrad_phi));
    phiphi = pagemtimes(reshape(phi,[],1,7),reshape(phi,1,[],7));
    
    h = sqrt(max(sum(d.^2,1)));
    Pe = mk*norm(problem.beta(center))*h/(2*problem.epsilon(center));
    tau = h/(2*norm(problem.beta(center)));
    if Pe < 1
        tau = mk*h^2/(4*problem.epsilon(center));
    end
    
    if SUPG == false
        tau = 0;
    end
    
    contrib_diff_tot = 2*area*squeeze(pagemtimes(omega .* epsilonpts, permute(gradgrad,[3,1,2])));
    contrib_conv_tot = 2*area*(omega.*betaTgrad) * phi';

    contrib_reaz_tot = 2*area*squeeze(pagemtimes(omega.*sigmapts, permute(phiphi,[3,1,2])));
    if MassLumping
        contrib_reaz_tot = diag(sum(contrib_reaz_tot,1));
    end
    
    contrib_massa_tot = 2*area*squeeze(pagemtimes(omega.*problem.rho(pts), permute(phiphi,[3,1,2])));

    local_laplacian = sum((invB*invB').*hessian_phi,[1,2]);
    
    for j=1:length(dof)
        jj = geom.pivot.pivot(dof(j));
        if jj > 0
            for k = 1:length(dof)
                contrib_diff = contrib_diff_tot(k,j);
                contrib_conv = contrib_conv_tot(k,j);
                contrib_reaz = contrib_reaz_tot(k,j);
                contrib_massa = contrib_massa_tot(k,j);

                contrib_SUPG = 2*area*tau*omega * (betaTgrad(k,:) .* betaTgrad(j,:))';
                contrib_SUPG = contrib_SUPG + 2*area*tau*omega*(local_laplacian(k) * betaTgrad(j,:).*epsilonpts)';
                kk = geom.pivot.pivot(dof(k));
                if kk > 0
                    A(jj,kk) = A(jj,kk) + contrib_diff + contrib_conv + contrib_reaz + contrib_SUPG;
                    mat_massa(jj,kk) = mat_massa(jj,kk) + contrib_massa;
                else
                    A_dirichlet(jj,-kk) = A_dirichlet(jj,-kk) + contrib_diff + contrib_conv + contrib_reaz + contrib_SUPG;
                    mat_massa_dirichlet(jj,-kk) = mat_massa_dirichlet(jj,-kk) + contrib_massa;
                end
            end
            b(jj) = b(jj) + ...
                2*area*(omega.*fpts)*phi(j,:)' + ...
                2*area*tau*(omega.*fpts)*betaTgrad(j,:)';
        end
    end
end

for j = 1:length(geom.pivot.Di)
    j = geom.pivot.Di(j,:);
    jj = -geom.pivot.pivot(j(1));
    coords = geom.elements.coordinates(j(1),:)';
    u_dirichlet(jj) = problem.bordo_dirichlet(coords,j(2)); %j(2) is marker
end

dof_border_position = logical(1:size(geom.elements.borders,2));
dof_border_position([3,4]) = 0;
for e = 1:length(geom.pivot.Ne)
    e = geom.pivot.Ne(e,:);
    points_idx = geom.elements.borders(e(1),dof_border_position); %[start node, end node, middle node]
    coords = geom.elements.coordinates(points_idx,:)';
    
    border_length = norm(coords(:,1)-coords(:,2));
    
    %from csi_1D to actual points along border and then get neumann values
    pts = coords(:,1) + (coords(:,2)-coords(:,1))*csi_1D;
    gamma_Ne = problem.bordo_neumann(pts, e(2)); %e(2) is marker
    
    % Numerical integration of the three basis functions
    contribs = border_length * (omega_1D .* gamma_Ne) * P_1D';
    
    % Accumulate
    jj = geom.pivot.pivot(points_idx);
    b_neumann(jj(jj > 0)) = b_neumann(jj(jj > 0)) + contribs(jj > 0)';
end


end