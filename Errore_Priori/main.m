clc
close all
clear all

global geom

epsilon = @(x) 1;
betax = @(x) 0;
betay = @(x) 0;
gamma = @(x) 0;
f = @(P) 32*(P(1)*(1-P(1))+P(2)*(1-P(2)));
bordo_dirichlet = @(P, marker) 0;
bordo_neumann = @(P, marker) -16*(P(:,2)-P(:,2).^2);

true_sol_handle = @(P) 0+16*P(1,:,:).*(1-P(1,:,:)).*P(2,:,:).*(1-P(2,:,:));



area_ax = linspace(0.005, 0.0001);
h_ax = 0*area_ax;
errors_ax = 0*area_ax;
tic
for i = 1:length(area_ax)
    generate_triangulation(area_ax(i));
    
    edge_start = geom.elements.coordinates(geom.elements.borders(:,1),:);
    edge_end = geom.elements.coordinates(geom.elements.borders(:,2),:);
    max_edge_length = sqrt(max(sum((edge_start-edge_end).^2,2)));
    h_ax(i) = max_edge_length;
    
    [A,b,A_dirichlet,u_dirichlet,b_neumann] = assembla(epsilon,betax,betay,gamma,f, bordo_dirichlet, bordo_neumann);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    
    [pippo] = calcola_errore_priori(true_sol_handle,u);
    errors_ax(i) = pippo;
    
%     sol = zeros(length(geom.elements.coordinates),1);
%     sol(geom.pivot.pivot > 0) = u;
%     sol(geom.pivot.pivot < 0) = u_dirichlet;
%     X = geom.elements.coordinates(:,1);
%     Y = geom.elements.coordinates(:,2);
%     true_sol = true_sol_handle([X,Y]')';
%     errors_ax(i) = norm(sol-true_sol,Inf);

end
toc
lm = fitlm(log(h_ax), log(errors_ax))
lm.coefCI
lm.plot