function [h_ax,errors_ax] = convergenzaErrorePriori(...
    Pk,SUPG,MassLumping,...
    true_sol_handle,grad_true_sol_handle,...
    nomeFile...
    )
global geom

area_ax = logspace(log10(0.1), log10(0.0001));
area_ax = area_ax(5:end);
h_ax = sqrt(area_ax);
errors_ax = zeros(length(area_ax),3);
for i = 1:length(area_ax)
    filename = append('Triangolazioni/',nomeFile,'/',num2str(area_ax(i)),'.mat');
    load(filename,'geom');
    if strcmp(Pk,'P2')
        prepP2()
    end
    
    [A,b,A_dirichlet,u_dirichlet,b_neumann,~,~] = assemblaEllittico(Pk,SUPG,MassLumping,false);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    
    sol = zeros(length(geom.elements.coordinates),1);
    sol(geom.pivot.pivot > 0) = u;
    sol(geom.pivot.pivot < 0) = u_dirichlet;
    
    [err_inf, err_L2, err_H1] = calcolaErrorePriori(true_sol_handle,grad_true_sol_handle,sol,Pk);
    errors_ax(i,:) = [err_inf, err_L2, err_H1];
end
end

