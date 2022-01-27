clear all
close all
clc
if(~exist('assemblaEllittico'))
    addpath('Funzioni')
end

global problem
global manifold
global geom

problem.epsilon = @(x) x(1,:)*0+1;
problem.beta = @(x) x*0;
problem.sigma = @(x) x(1,:)*0;
problem.f = @(x,t) x(1,:)*0;
problem.rho = @(x) x(1,:)*0+1;

manifold.name = 'sphere';
manifold.parameters.r = 1;

manifold_geodiff();


N = 3;
eig_true = [];
for k = 0:N
    for m = 1:2*k+1
        eig_true(end+1) = k*(k+1);
    end
end
eig_true
%%
area_ax = logspace(log10(0.01), log10(0.0001));
%area_ax = area_ax(10:30);
h_ax = sqrt(area_ax);

eig_approx = zeros(length(h_ax),length(eig_true));
for m = 1:length(h_ax)
    [m,round(1/h_ax(m))]
    manifold_triangulator(round(1/h_ax(m)));
    [A,~,M] = manifold_assemblaEllittico('P1',true);    
    eig_approx(m,:) = sort(eigs(A,M,length(eig_true),'smallestabs'))';
end

%% Eigenvectors
[V,D] = eigs(A,M,length(eig_true),'smallestabs');
Vtilde = V(geom.pivot.pivot,:);
map_surf = manifold.numeric.P(geom.elements.coordinates');
figure(1)
clf
tiledlayout(4,4, 'Padding', 'none', 'TileSpacing', 'none'); 
for i = 1:length(eig_true)
    nexttile
    trisurf(geom.elements.triangles, map_surf(1,:), map_surf(2,:), map_surf(3,:), Vtilde(:,i))
    title("\lambda = "+num2str(eig_true(i)))
    shading interp
    caxis([min(Vtilde,[],'all'),max(Vtilde,[],'all')])
    axis equal
    axis off
end


%% Export to LATEX
%writematrix([h_ax', errors],'convergenza_spettro_sfera.csv')

%%
close all
clc
errors = abs(eig_true-eig_approx);
p = loglog(h_ax, errors);
for l = unique(eig_true)
    if l > 0
        set(p(eig_true == l),'Color', rand(1,3));
    end
end
legend(cellstr(num2str(eig_true', 'eig =%-d')),'Location','southeast')
delete(p(1))