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
grad_true_sol_handle = @(P) [   16*P(2,:,:).*(1-P(2,:,:)).*(1-2*P(1,:,:)); 
                                16*P(1,:,:).*(1-P(1,:,:)).*(1-2*P(2,:,:))];


area_ax = logspace(log10(0.08), log10(0.0002), 12);
h_ax = sqrt(area_ax);
errors_ax = 0*area_ax;
tic
for i = 1:length(area_ax)
    generate_triangulation(area_ax(i));
        
    [A,b,A_dirichlet,u_dirichlet,b_neumann] = assembla(epsilon,betax,betay,gamma,f, bordo_dirichlet, bordo_neumann);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    
    sol = zeros(length(geom.elements.coordinates),1);
    sol(geom.pivot.pivot > 0) = u;
    sol(geom.pivot.pivot < 0) = u_dirichlet;
    
    
    [tot_err_inf, tot_err_L2, tot_err_H1] = calcola_errore_priori(true_sol_handle,grad_true_sol_handle,sol);
    errors_ax(i) = tot_err_inf;

end
toc
p = polyfit(log(h_ax), log(errors_ax),1)
f = polyval(p,log(h_ax));
plot(log(h_ax),log(errors_ax),'o',log(h_ax),f,'-')
legend('data','linear fit')
%lm = fitlm(log(h_ax), log(errors_ax))
%lm.coefCI
%lm.plot