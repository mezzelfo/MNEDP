clear all
close all
clc
if(~exist('assemblaEllittico'))
     addpath('Funzioni')
end

%% %TODO: Chiedere a Berrone se bisogna tenere u_true fissa o f fissa.
global problem
true_sol_handle = @(x) 5+16*x(1,:).*(1-x(1,:)).*x(2,:).*(1-x(2,:));
grad_true_sol_handle = @(x) [   16*x(2,:).*(1-x(2,:)).*(1-2*x(1,:));
    16*x(1,:).*(1-x(1,:)).*(1-2*x(2,:))];

problem.beta = @(x) x*0+[1;1];
problem.sigma = @(x) x(1,:)*0;
problem.bordo_dirichlet = @(x, marker) 5;

epsilon_ax = logspace(1,-6,10);
noSUPG_res = [];
yesSUPG_res = [];

for epsilon = epsilon_ax
    tic
    problem.epsilon = @(x) x(1,:)*0+epsilon;
    problem.f = @(x) epsilon*32*(x(1,:).*(1-x(1,:))+x(2,:).*(1-x(2,:))) + dot(problem.beta(x), grad_true_sol_handle(x));
    problem.bordo_neumann = @(x, marker) epsilon*(-16*(x(2,:)-x(2,:).^2));

    [h_ax,errors_ax_P1_noSUPG] = convergenzaErrorePriori(...
        'P1',false,false,...
        true_sol_handle,grad_true_sol_handle,...
        'QuadratoMisto'...
        );

    noSUPG_res(end+1,:,:) = errors_ax_P1_noSUPG;

    [h_ax,errors_ax_P1_yesSUPG] = convergenzaErrorePriori(...
        'P1',true,false,...
        true_sol_handle,grad_true_sol_handle,...
        'QuadratoMisto'...
        );

    yesSUPG_res(end+1,:,:) = errors_ax_P1_yesSUPG;
    toc
end

polyfit(log(h_ax), log(yesSUPG_res(end,:,2)),1)

%% Plot

figure(1)
subplot(1,2,1)
loglog(h_ax,noSUPG_res(:,:,2)) 
title('senza SUPG')
xlabel('h')
ylabel('Errore L2')

subplot(1,2,2)
loglog(h_ax,yesSUPG_res(:,:,2))
title('con SUPG')
xlabel('h')
ylabel('Errore L2')


sgtitle('Errore L2 FEM P1: SUPG vs no SUPG')

figure(2)
subplot(1,2,1)
loglog(h_ax,noSUPG_res(:,:,3)) 
title('senza SUPG')
xlabel('h')
ylabel('Errore H1')

subplot(1,2,2)
loglog(h_ax,yesSUPG_res(:,:,3))
title('con SUPG')
xlabel('h')
ylabel('Errore H1')


sgtitle('Errore H1 FEM P1: SUPG vs no SUPG')