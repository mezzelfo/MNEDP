clc
close all
clear all

true_sol_handle = @(P) 5+16*P(1,:).*(1-P(1,:)).*P(2,:).*(1-P(2,:));
grad_true_sol_handle = @(P) [   16*P(2,:).*(1-P(2,:)).*(1-2*P(1,:));
    16*P(1,:).*(1-P(1,:)).*(1-2*P(2,:))];

M = 0.0001;
epsilon = @(x) x(1,:)*0+M;
beta = @(x) x*0+M*[100000;1];
f = @(P) M*32*(P(1,:).*(1-P(1,:))+P(2,:).*(1-P(2,:))) + dot(beta(P), grad_true_sol_handle(P));
bordo_dirichlet = @(P, marker) 5;
bordo_neumann = @(P, marker) M*(-16*(P(2,:)-P(2,:).^2));

[h_ax,errors_ax_P1_yesSUPG] = convergenzaErrorePriori(...
    epsilon,beta,f, bordo_dirichlet, bordo_neumann,...
    'P1',true,...
    true_sol_handle,grad_true_sol_handle...
    );

[h_ax,errors_ax_P1_noSUPG] = convergenzaErrorePriori(...
    epsilon,beta,f, bordo_dirichlet, bordo_neumann,...
    'P1',false,...
    true_sol_handle,grad_true_sol_handle...
    );

loglog(h_ax,errors_ax_P1_yesSUPG(:,2),'-o',h_ax,errors_ax_P1_noSUPG(:,2),'-o')
legend('with SUPG','without SUPG');


% beta_ax = logspace(-1,6,10*2);
% yesSUPG_ax = [];
% noSUPG_ax = [];
%
% for b = beta_ax
%     tic
%     epsilon = @(x) x(1,:)*0+1;
%     beta = @(x) x*0+[b;1];
%     f = @(P) 32*(P(1,:).*(1-P(1,:))+P(2,:).*(1-P(2,:))) + dot(beta(P), grad_true_sol_handle(P));
%     bordo_dirichlet = @(P, marker) 5;
%     bordo_neumann = @(P, marker) -16*(P(2,:)-P(2,:).^2);
%
%     [h_ax,errors_ax_P1_yesSUPG] = convergenzaErrorePriori(...
%         epsilon,beta,f, bordo_dirichlet, bordo_neumann,...
%         'P2',true,...
%         true_sol_handle,grad_true_sol_handle...
%         );
%
%     p = polyfit(log(h_ax), log(errors_ax_P1_yesSUPG(:,2)),1);
%     yesSUPG_ax(end+1) = p(1);
%
%     [h_ax,errors_ax_P1_noSUPG] = convergenzaErrorePriori(...
%         epsilon,beta,f, bordo_dirichlet, bordo_neumann,...
%         'P2',false,...
%         true_sol_handle,grad_true_sol_handle...
%         );
%
%     p = polyfit(log(h_ax), log(errors_ax_P1_noSUPG(:,2)),1);
%     noSUPG_ax(end+1) = p(1);
%     toc
% end
%
% semilogx(beta_ax,yesSUPG_ax,'-o',beta_ax,noSUPG_ax,'-o')
% legend('with SUPG','without SUPG')
% xlabel('first component of beta')
% ylabel('Estimated asymptotic order of convergence')
