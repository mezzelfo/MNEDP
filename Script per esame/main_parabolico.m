clear all
close all
clc
if(~exist('assemblaEllittico'))
     addpath('Funzioni')
end

global geom
global problem

generaTriangolazione(0.006, [-1,-1;1,-1;1,1;-1,1], [3,3,3,3])
%%

problem.epsilon = @(x) x(1,:)*0+1;
problem.beta = @(x) [0*x(1,:);0*x(2,:)];
problem.sigma = @(x) 0*x(1,:);
problem.f = @(x,t) 0*x(1,:);
problem.bordo_dirichlet = @(x,t, marker) 0;
problem.bordo_neumann = @(x,t, marker) 0;
problem.rho = @(x) 0*x(1,:) + 1;
problem.iniziale = @(x) exp(-vecnorm([x(1,:);x(2,:)]).^2);

n_steps = 100;
T = 1;
Pk = 'P2';

if isequal(Pk,'P2') && (geom.nelements.nVertexes == length(geom.elements.coordinates)) 
    prepP2
end
[u,uD] = assemblaParabolico(Pk,T,n_steps);

% Assembla la soluzione
utilde = zeros(geom.nelements.nVertexes,n_steps);

for t = 1:n_steps
    for v = 1:geom.nelements.nVertexes
        ii = geom.pivot.pivot(v);
        if ii > 0 % Dof
            utilde(v,t) = u(ii,t);
        else
            utilde(v,t) = uD(-ii,t);
        end
    end
end

%%

nV = geom.nelements.nVertexes;
for t = 1:n_steps
    trisurf(geom.elements.triangles(:,1:3), geom.elements.coordinates(1:nV,1), geom.elements.coordinates(1:nV,2), utilde(:,t))
    shading interp
    axis equal
    axis([-1 1 -1 1 -2 2])
    pause(0.03)
    if t ~= n_steps
        clf
    end
end
