clc
close all
clear all

global geom

true_sol_handle = @(P) 0+16*P(1,:).*(1-P(1,:)).*P(2,:).*(1-P(2,:));
grad_true_sol_handle = @(P) [   16*P(2,:).*(1-P(2,:)).*(1-2*P(1,:)); 
                                16*P(1,:).*(1-P(1,:)).*(1-2*P(2,:))];

epsilon = @(x) 1;
beta = @(x) [10000;2];
f = @(P) 32*(P(1,:)*(1-P(1,:))+P(2,:)*(1-P(2,:))) + beta(P)'*grad_true_sol_handle(P);
bordo_dirichlet = @(P, marker) 0;
bordo_neumann = @(P, marker) -16*(P(2,:)-P(2,:).^2);


area_ax = logspace(log10(0.08), log10(0.0002), 12);
h_ax = sqrt(area_ax);
errors_ax = 0*area_ax;
tic
for i = 1:length(area_ax)
    
    filename = append('triangolazione_',num2str(area_ax(i)),'.mat');
    %generate_triangulation(area_ax(i));
    %save(filename,'geom')
    load(filename);
        
    [A,b,A_dirichlet,u_dirichlet,b_neumann] = assembla(epsilon,beta,f, bordo_dirichlet, bordo_neumann);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    
    sol = zeros(length(geom.elements.coordinates),1);
    sol(geom.pivot.pivot > 0) = u;
    sol(geom.pivot.pivot < 0) = u_dirichlet;
    
    
    [tot_err_inf, tot_err_L2, tot_err_H1] = calcola_errore_priori(true_sol_handle,grad_true_sol_handle,sol);
    errors_ax(i) = tot_err_L2;

end
toc
p = polyfit(log(h_ax), log(errors_ax),1)
f = polyval(p,log(h_ax));
plot(log(h_ax),log(errors_ax),'o',log(h_ax),f,'-')
legend('data','linear fit')
%lm = fitlm(log(h_ax), log(errors_ax))
%lm.coefCI
%lm.plot