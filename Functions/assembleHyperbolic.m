function [u,uD] = assembleHyperbolic(T,n_steps,beta_Newmark,fzero)

% Solves an hyperbolic diffusion-convection-reaction problem on a 2-manifold
% embedded in R^3 with polynomial finite elements and Newmark advancement method.
% The polynomial order is given by "Pk", the time horizon by "T", the number
% of steps "n_steps" and the Newmark method's parameter by "beta_Newmark".
% Given the structs "geom" built with "assembleManifold", "geom" with
% "triangulateChart" and "problem" containing
% - rho: mass density function (constant in time)
% - epsilon: diffusivity function (constant in time)
% - beta: convective velocity field (constant in time)
% - sigma: reaction function (constant in time)
% - f: function describing sources and sinks (may be time dependent)
% - boundary_D: function describing dirichlet boundary conditions (constant in time)
% - initial: starting solution
% - initial_v: starting velocity
% returns the sequence of solutions in the interior points "u" and on the
% Dirichlet boundary "uD". If "f" is zero then the option "fzero"=true will
% greatly reduce the computational cost.

global geom;
global manifold;
global problem;

% Obtain the stiffness and mass matrices together with the vector b at time 0
f_time = problem.f;
problem.f = @(x) f_time(x,0);
[A,AD,b_old_old,uD,M,~] = assembleElliptic(true);
problem.f = f_time;

% Time step length
delta_t = T/n_steps;

b_old = zeros(size(b_old_old));
b_new = zeros(size(b_old_old));

% Preallocate the solution
u = zeros(length(b_new),n_steps+1);

% Quadrature nodes in the reference triangle
[zita,csi,eta,omega] = int_nodes_weights();

% Basis functions in the reference triangle
phi = [csi; eta; zita];

% Interpolate the initial condition
for n = 1:geom.nelements.nVertexes
    nn = geom.pivot.pivot(n);
    if nn > 0
        u(nn,1) = problem.initial(geom.elements.coordinates(n,:)');
    end
end
% Interpolate the initial velocity
v = zeros(length(b_new),1);
for n = 1:geom.nelements.nVertexes
    nn = geom.pivot.pivot(n);
    if nn > 0
        v(nn) = problem.initial_v(geom.elements.coordinates(n,:)');
    end
end

% Compute the first step with explicit Newton
u(:,2) = u(:,1) + delta_t*v;

% Assemble b_old
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
    b_old(ii(ii > 0)) = b_old(ii(ii > 0)) + 2*area*sum(omega.*sqrtdetg.*f.*phi(ii > 0,:),2);
end

% Precompute useful matrices
m_Newmark = M + beta_Newmark*delta_t^2*A;
m_Newmark2 = 2*M - (1-2*beta_Newmark)*delta_t^2*A;
ADuD = AD*uD;

% Time advancement with Newmark's method
for t = 1:n_steps-1
    if fzero % Computational speed
        u(:,t+2) = m_Newmark\(m_Newmark2*u(:,t+1) - m_Newmark*u(:,t) - delta_t^2*ADuD);
    else
        for e = 1:geom.nelements.nTriangles
            indices = geom.elements.triangles(e,:); % Vertices indices
            coords = geom.elements.coordinates(indices,:)'; % Vertices coordinates
            d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
            area = geom.support.TInfo(e).Area;
            B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
            pts = coords(:,end)+B*[csi; eta];
            
            f = problem.f(pts,t*delta_t);
            sqrtdetg = manifold.numeric.sqrtdetg(pts);
            ii = geom.pivot.pivot(indices);
            b_new(ii(ii > 0)) = b_new(ii(ii > 0)) + 2*area*sum(omega.*sqrtdetg.*f.*phi(ii > 0,:),2);
        end
        u(:,t+2) = m_Newmark\(m_Newmark2*u(:,t+1) - m_Newmark*u(:,t) - delta_t^2*ADuD + ...
            delta_t^2*(beta_Newmark*b_new + (1-2*beta_Newmark)*b_old + beta_Newmark*b_old_old));
        b_old_old = b_old;
        b_old = b_new;
        b_new = 0*b_new;
    end
end
end