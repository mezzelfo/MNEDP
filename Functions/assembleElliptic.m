function [A,AD,b,uD,M,MD] = assembleElliptic(mass)
% Solves an elliptic diffusion-convection-reaction problem on a 2-manifold
% embedded in R^3 with polynomial finite elements.
% Given the structs "geom" built with "assembleManifold", "geom" with
% "triangulateChart" and "problem" containing
% - rho: mass density (needed only if mass = true)
% - epsilon: diffusivity function
% - beta: convective velocity field
% - sigma: reaction function
% - f: function describing sources and sinks
% - boundary_D: function describing dirichlet boundary conditions
% Returns the stiffness matrices "A", "AD" (for dirichlet points), the
% vectors "b" such that the Galerkin homogeneous solution "u0" is a solution of the
% linear system A*u0 = b - AD*uD and the solution on the dirichlet boundary
% "uD".
% If mass = true the function return also the mass matrices M, MD

global geom;
global manifold;
global problem;

ndof = max(geom.pivot.pivot); % Number of degrees of freedom

A = spalloc(ndof,ndof,10*ndof); % Preallocate the stiffnes matrix
if mass
    M = spalloc(ndof,ndof,10*ndof);
    if any(manifold.dirichlet_borders)
        MD = spalloc(ndof,-min(geom.pivot.pivot),-10*min(geom.pivot.pivot));
    else
        MD = 0;
    end
else
    M = NaN;
    MD = NaN;
end
b = zeros(ndof,1);
% Preallocate the Dirichlet stiffness matrix
if any(manifold.dirichlet_borders)
    AD = spalloc(ndof, -min(geom.pivot.pivot),-10*min(geom.pivot.pivot));
    uD = zeros(-min(geom.pivot.pivot),1);
else
    AD = 0;
    uD = 0;
end

% Nodes and weights for quadrature inside the reference triangle
[zita,csi,eta,omega] = int_nodes_weights();
n_quad = length(omega);

% Basis functions and gradients inside the reference triangle
phi = [csi; eta; zita];
grad_phi = [[1 0];[0 1];[-1 -1]]';
grad_phi = repmat(grad_phi,1,1,length(omega));

Ginv = zeros(2,2,length(csi));

for e = 1:geom.nelements.nTriangles  % Cycle over the triangles
    points_idx = geom.elements.triangles(e,:);  % Vertices indices
    coords = geom.elements.coordinates(points_idx,:)'; % Vertices coordinates
    d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
    area = geom.support.TInfo(e).Area; % Triangle area
    B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
    invB = inv(B);
    pts = coords(:,end)+B*[csi; eta]; % Quadrature points coordinates
    sqrtdetg = manifold.numeric.sqrtdetg(pts); % Volume form in the quadrature points
    
    % Inverse metric in the quadrature points
    for c = 1:n_quad
        Ginv(:,:,c) = manifold.numeric.ginv(pts(:,c));
    end
    epsilon = problem.epsilon(pts); % Diffusivity
    beta = problem.beta(pts); % Velocity field
    sigma = problem.sigma(pts); % Reaction term
    f = problem.f(pts); % Forcing term
    
    BinvTgrad_phi = pagemtimes(invB,'transpose',grad_phi,'none');
    Ginvgrad = pagemtimes(Ginv,BinvTgrad_phi);
    gradGinvgrad = pagemtimes(BinvTgrad_phi,'transpose',Ginvgrad,'none');
    
    betaTgrad = squeeze(sum(reshape(beta,[2,1,7]).*BinvTgrad_phi));
    phiphi = pagemtimes(reshape(phi,[],1,7),reshape(phi,1,[],7));
    
    % Diffusion, convection and reaction contributions to the stiffness matrix
    diffusion = 2*area*squeeze(pagemtimes(omega.*sqrtdetg.*epsilon, permute(gradGinvgrad,[3,1,2])));
    convection = 2*area*(omega.*sqrtdetg.*betaTgrad) * phi';
    reaction = 2*area*squeeze(pagemtimes(omega.*sqrtdetg.*sigma, permute(phiphi,[3,1,2])));
    total = diffusion+convection+reaction;
    
    if mass % Compute the triangle's contribution to the mass matrix
        rho = problem.rho(pts);
        mass_contribution = 2*area*squeeze(pagemtimes(omega.*sqrtdetg.*rho, permute(phiphi,[3,1,2])));
    end
    
    for j=1:size(phi,1)
        jj = geom.pivot.pivot(points_idx(j));
        if jj > 0
            for k = 1:size(phi,1)
                kk = geom.pivot.pivot(points_idx(k));
                if kk > 0
                    A(jj,kk) = A(jj,kk) + total(k,j);
                    if mass
                        M(jj,kk) = M(jj,kk) + mass_contribution(k,j);
                    end
                else % Dirichlet
                    AD(jj,-kk) = AD(jj,-kk) + total(k,j);
                    if mass
                        MD(jj,-kk) = MD(jj,-kk) + mass_contribution(k,j);
                    end
                end
            end
            b(jj) = b(jj) + 2*area*(omega.*sqrtdetg.*f)*phi(j,:)'; % Termine noto
        else % Dirichlet
            uD(-jj) = problem.boundary_D(coords(:,j));
        end
    end
end
end