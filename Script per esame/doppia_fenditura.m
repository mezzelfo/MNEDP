clear all
close all
clc
if(~exist('assemblaEllittico'))
    addpath('Funzioni')
end

global problem
global geom
global manifold


manifold.uv_bounds = [-15,15;-12,12];
manifold.dirichlet_borders = [0,0,0,0];
manifold.numeric.P = @(x) x;
N = 500;
manifold_triangulator(N)

geom.pivot.Di = [];
geom.pivot.Ne = [];
geom.elements.borders = [];

%% Dati del problema
problem.epsilon = @(x) x(1,:)*0-1;
problem.rho = @(x) x(1,:)*0+1i;
problem.beta = @(x) x*0;
problem.sigma = @(x) 10000*(abs(x(1,:)) <= 0.25).*(-1+(abs(x(2,:)-0.5) <= 0.25)+(abs(x(2,:)+0.5) <= 0.25));
problem.f = @(x,t) x(1,:)*0;
x0 = -7;
y0 = 0;
delta_x = 0.7;
delta_y = 1.5;
kx0 = 20;
ky0 = 0;
problem.iniziale = @(x) 1/(2*delta_x^2*pi)^(1/4)*1/(2*delta_y^2*pi)^(1/4)*exp(-((x(1,:)-x0)./(2*delta_x)).^2).*exp(-((x(2,:)-y0)./(2*delta_y)).^2).*exp(1i * (kx0*x(1,:) + ky0*x(2,:)));
problem.bordo_dirichlet = @(x,t,marker) x(1,:)*0;
problem.bordo_neumann = @(x,t,marker) x(1,:)*0;

%% Risoluzione
n_steps = 200;
T = 1;
tic
[u,uD] = assemblaParabolico('P1',T,n_steps,false,true,true);
toc
% Assembla la soluzione
utilde = zeros(geom.nelements.nVertexes,n_steps+1);
utilde(geom.pivot.pivot > 0,:) = u;
utilde(geom.pivot.pivot < 0,:) = uD;

%% SAVE utilde for BLENDER processing
%save('doppia_fenditura_utilde.mat','utilde');
for t = 1:n_steps+1
    img = abs(reshape(utilde(:,t),N,N)).^2;
    imwrite(img,"Fenditura_OUTPUT/img_"+num2str(t)+".png")
end

%% Plot soluzione
figure(1)
for t = 1:n_steps+1
    clf
    trisurf(...
        geom.elements.triangles,...
        geom.elements.coordinates(:,1),...
        geom.elements.coordinates(:,2),...
        abs(utilde(:,t)).^2)
    view(0,90)
    shading interp
    title(num2str(t*(T/n_steps)))
    pause(0.001)
end

%% Plot potenziale
% figure(1)
% trisurf(...
%     geom.elements.triangles,...
%     geom.elements.coordinates(:,1),...
%     geom.elements.coordinates(:,2),...
%     problem.sigma(geom.elements.coordinates'))
% axis square
% view(0,90)

%% Conservazione della probabilità
% integrale = [];
% for t = 1:n_steps+1
%     [~, tot_err_L2, ~] = calcolaErrorePriori(@(x) x(1,:)*0, @(x) x*0, abs(utilde(:,t)), 'P1');
%     integrale(end+1) = tot_err_L2;
% end
% plot(integrale)
% ylim([0,1])

%% Plot Soluzione
% borders = geom.pivot.Ne(geom.pivot.Ne(:,2) == 100,1);
% pts_idx = geom.elements.borders(borders,1:2);
% pts_idx = unique(pts_idx(:));
% pts = geom.elements.coordinates(pts_idx,2);
% [pts,idx] = sort(pts);
% pts_idx = pts_idx(idx);
%
%
% nv = geom.nelements.nVertexes;
% M = max(max(abs(utilde(pts_idx,:))));

%%
figure(1)
clf
pts_mask = abs(geom.elements.coordinates(:,1) - 7.0040) < 0.001;
pts = geom.elements.coordinates(pts_mask, 2);
[pts,idx] = sort(pts);
polipo = utilde(pts_mask,:);
t_ax = 130:160;
param = linspace(0,1,length(t_ax));
C = ([0,0,1]'*param + [1,0,0]'*(1-param))';
axes('ColorOrder',C,'NextPlot','replacechildren')
plot(pts,abs(polipo(idx,t_ax)).^2,'LineWidth',1)
colormap(C)
caxis([min(t_ax),max(t_ax)])
c = colorbar;
c.Label.String = 'passo temporale';
ylim([0,0.035])
xlabel('y')
ylabel('|\Psi(6,y)|^2')



