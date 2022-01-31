clear all
close all
clc
if(~exist('assembleElliptic'))
    addpath('Functions')
end


global geom;
global manifold;
global problem;


manifold.name = 'sphere';
manifold.parameters.r = 1;

assembleManifold();
triangulateChart(40,40);
%%

problem.epsilon = @(x) 0*x(1,:)+1;
problem.beta = @(x) 0*x;
problem.sigma = @(x) 0*x(1,:);
problem.f = @(x) sin(3*x(1,:)).*cos(10*x(2,:));
problem.boundary_D = @(x) sin(2*x(1,:));


[A,AD,b,uD,~,~] = assembleElliptic(false);
% To grant uniqueness to the solution one might need other conditions
% A = [A;ones(1,size(A,2))];
% b = [b;0];

u0 = A\(b-AD*uD);
usol = zeros(geom.nelements.nVertexes,1);
for i = 1:geom.nelements.nVertexes
    ii = geom.pivot.pivot(i);
    if ii > 0
        usol(i) = u0(ii);
    else
        usol(i) = uD(-ii);
    end
end

%% Plot solution
map_surf = manifold.numeric.P(geom.elements.coordinates');
subplot(1,2,1)
trisurf(geom.elements.triangles, map_surf(1,:), map_surf(2,:), map_surf(3,:), usol)
axis equal
shading interp
title("Solution on the manifold")
subplot(1,2,2)
trisurf(geom.elements.triangles, geom.elements.coordinates(:,1),...
    geom.elements.coordinates(:,2), usol)
axis equal
shading interp
title("Solution in the chart")