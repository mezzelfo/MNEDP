function [h_ax,errors_ax] = convergenzaErrorePriori(...
    epsilon,beta,f, bordo_dirichlet, bordo_neumann,...
    Pk,SUPG,...
    true_sol_handle,grad_true_sol_handle...
    )
global geom
area_ax = logspace(log10(0.1), log10(0.0001));
area_ax = area_ax(1:10:end);
h_ax = sqrt(area_ax);
errors_ax = zeros(length(area_ax),3);
for i = 1:length(area_ax)
    
    filename = append('../Triangolazioni/triangolazione_',num2str(area_ax(i)),'.mat');
    %generate_triangulation(area_ax(i));
    %save(filename,'geom')
    load(filename,'geom');
    if strcmp(Pk,'P2')
        prepP2()
    end
    
    [A,b,A_dirichlet,u_dirichlet,b_neumann] = assembla(epsilon,beta,f,bordo_dirichlet,bordo_neumann,SUPG,Pk);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    
    sol = zeros(length(geom.elements.coordinates),1);
    sol(geom.pivot.pivot > 0) = u;
    sol(geom.pivot.pivot < 0) = u_dirichlet;
    
    [err_inf, err_L2, err_H1] = calcola_errore_priori(true_sol_handle,grad_true_sol_handle,sol,Pk);
    errors_ax(i,:) = [err_inf, err_L2, err_H1];
end
end

