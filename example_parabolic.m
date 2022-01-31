clear
close all
clc
if(~exist('assembleParabolic'))
    addpath('Functions')
end


global geom;
global manifold;
global problem;


manifold.name = 'sphere';
manifold.parameters.r = 1;
assembleManifold();
triangulateChart(50,50);
%%

problem.rho = @(x) 0*x(1,:) + 1;
problem.epsilon = @(x) 0*x(1,:)+1;
problem.beta = @(x) 0*x;
problem.sigma = @(x) 0*x(1,:);
problem.f = @(x,t) 0*x(2,:);
problem.boundary_D = @(x) 0*x(1,:);
problem.initial = @(x) exp(-x(2,:).^2/0.5);

T = 2;
n_steps = 60;
[u,uD] = assembleParabolic(T,n_steps,true);

usol = zeros(geom.nelements.nVertexes,n_steps+1);
for t = 1:n_steps+1
    for i = 1:geom.nelements.nVertexes
        ii = geom.pivot.pivot(i);
        if ii > 0
            usol(i,t) = u(ii,t);
        else
            usol(i,t) = uD(-ii);
        end
    end
end
%% Plot solution

map_surf = manifold.numeric.P(geom.elements.coordinates');
figure
pause
for t = 1:20
    pause(1)
    subplot(1,2,1)
    trisurf(geom.elements.triangles, map_surf(1,:), map_surf(2,:), map_surf(3,:), usol(:,t))
    axis equal
    shading interp
    title("Solution on the manifold")
    caxis([-0.2 0.5])
    subplot(1,2,2)
    trisurf(geom.elements.triangles, geom.elements.coordinates(:,1),...
        geom.elements.coordinates(:,2), usol(:,t))
    axis equal
    zlim([-0.5 1])
    caxis([-0.2 0.5])
    shading interp
    title("Solution in the chart")
    colorbar
end