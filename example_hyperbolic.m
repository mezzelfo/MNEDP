clear
close all
clc
if(~exist('assembleHyperbolic'))
    addpath('Functions')
end


global geom;
global manifold;
global problem;


manifold.name = 'torus';
manifold.parameters.r1 = 1;
manifold.parameters.r2 = 2;
assembleManifold();
triangulateChart(90,90);
%%

problem.rho = @(x) 0*x(1,:) + 1;
problem.epsilon = @(x) 0*x(1,:) + 1;
problem.beta = @(x) 0*x;
problem.sigma = @(x) 0*x(1,:);
problem.f = @(x,t) 0*x(2,:);
problem.boundary_D = @(x) 0*x(1,:);
problem.initial = @(x) exp(-vecnorm(x-[pi/2;pi/2]).^2*10) + exp(-vecnorm(x-[3*pi/2;3*pi/2]).^2*10);
problem.initial_v = @(x) 0*x(1,:);

T = 6;
n_steps = 600;
[u,uD] = assembleHyperbolic(T,n_steps,1/4,true);

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

gif = true;
map_surf = manifold.numeric.P(geom.elements.coordinates');
normal = zeros(size(map_surf));
for i = 1:geom.nelements.nVertexes
    normal(:,i) = 0.7*manifold.numeric.normal(geom.elements.coordinates(i,:)');
end
h = figure;
pause
for t = 1:3:n_steps+1
    pause(0.005)
    subplot(1,2,1)
    trisurf(geom.elements.triangles, map_surf(1,:) + usol(:,t)'.*normal(1,:),...
        map_surf(2,:) + usol(:,t)'.*normal(2,:), map_surf(3,:) + usol(:,t)'.*normal(3,:), usol(:,t))
    axis equal
    axis([-3.3 3.3 -3.3 3.3 -2 2])
    view(30,35)
    shading interp
    caxis([-0.1 0.3])
    title("Solution on the manifold")
    subplot(1,2,2)
    trisurf(geom.elements.triangles, geom.elements.coordinates(:,1),...
        geom.elements.coordinates(:,2), usol(:,t))
    axis equal
    axis([0 2*pi 0 2*pi -0.5 1])
    shading interp
    caxis([-0.1 0.3])
    title("Solution in the chart")
    colorbar
    
    if gif
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if t == 1
            imwrite(imind,cm,'torus_wave.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'torus_wave.gif','gif','WriteMode','append');
        end
    end
end

