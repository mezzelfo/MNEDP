function [u,uD] = assembleParabolic(T,n_steps,fzero)

% Solves a parabolic diffusion-convection-reaction problem on a 2-manifold
% embedded in R^3 with polynomial finite elements and Crank-Nicolson advancement method.
% The polynomial order is given by "Pk", the time horizon by "T" and the number
% of steps "n_steps"
% Given the structs "geom" built with "assembleManifold", "geom" with
% "triangulateChart" and "problem" containing
% - rho: mass density function (constant in time)
% - epsilon: diffusivity function (constant in time)
% - beta: convective velocity field (constant in time)
% - sigma: reaction function (constant in time)
% - f: function describing sources and sinks (may be time dependent)
% - boundary_D: function describing dirichlet boundary conditions (constant in time)
% - initial: starting solution
% returns the sequence of solutions in the interior points "u" and on the
% Dirichlet boundary "uD".
% If "f" is zero then the option "fzero"=true will greatly reduce the computational cost.

global geom;
global manifold;
global problem;

% Obtain the stiffness and mass matrices together with the vector b at time 0
f_time = problem.f;
problem.f = @(x) f_time(x,0); % Forcing function at time 0 needed for b_old
[A,AD,b_old,uD,M,~] = assembleElliptic(true);
problem.f = f_time;

delta_t = T/n_steps; % Time step length
b_new = zeros(length(b_old),1);

% Preallocate solution
u = zeros(length(b_new),n_steps+1);

% Quadrature nodes in the reference triangle
[zita,csi,eta,omega] = int_nodes_weights();

% Basis functions in the reference triangle
phi = [csi; eta; zita];

% Starting solution interpolation
for n = 1:geom.nelements.nVertexes
    nn = geom.pivot.pivot(n);
    if nn > 0
        u(nn,1) =  problem.initial(geom.elements.coordinates(n,:)');
    end
end

% Useful precomputations
Kp = M + delta_t/2*A; % M and A are constant in time
Km =  M - delta_t/2*A;
ADuD = AD*uD; % Both AD and uD are constant in time

% Time advancement with Crank-Nicolson method
for t = 1:n_steps
    % Compute the new forcing vector b_new
    b_new = 0*b_new;
    if ~fzero
        for e = 1:geom.nelements.nTriangles
            indices = geom.elements.triangles(e,:); % Vertices indices
            coords = geom.elements.coordinates(indices,:)'; % Vertices coordinates
            d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
            area = geom.support.TInfo(e).Area; % Triangle area
            B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
            pts = coords(:,end)+B*[csi; eta]; % Quadrature points
            
            f = problem.f(pts,t*delta_t);
            sqrtdetg = manifold.numeric.sqrtdetg(pts);
            ii = geom.pivot.pivot(indices);
            b_new(ii(ii > 0)) = b_new(ii(ii > 0)) + 2*area*(sqrtdetg.*f.*phi(ii > 0,:))*omega';
            
        end
    end
    % Compute the solution in the next time
    u(:,t+1) = Kp\(Km*u(:,t)+delta_t/2*(b_new + b_old - 2*ADuD));
    b_old = b_new;
end
end