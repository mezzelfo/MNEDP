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

problem.beta = @(x) x*0+[0;0];
problem.sigma = @(x) x(1,:)*0 + 1;
problem.bordo_dirichlet = @(P, marker) 5;

epsilon_ax = logspace(1,-6,10);
noML_res = [];
yesML_res = [];

for epsilon = epsilon_ax
    tic
    problem.epsilon = @(x) x(1,:)*0+epsilon;
    problem.f = @(x) epsilon*32*(x(1,:).*(1-x(1,:))+x(2,:).*(1-x(2,:))) + problem.sigma(x).*true_sol_handle(x);
    problem.bordo_neumann = @(x, marker) epsilon*(-16*(x(2,:)-x(2,:).^2));
    
    %%%Plot... non riusciamo a farla oscillare
    global geom
    generaTriangolazione(0.003, [0,0;1,0;1,1;0,1], [3,3,3,3])
    [A,b,A_dirichlet,u_dirichlet,b_neumann,s,se] = assemblaEllittico('P1',false,true,false);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    utilde = zeros(geom.nelements.nVertexes,1);
    utilde(geom.pivot.pivot>0) = u;
    utilde(geom.pivot.pivot<0) = u_dirichlet;
    trisurf(geom.elements.triangles,geom.elements.coordinates(:,1),geom.elements.coordinates(:,2),utilde)
    return
    %%%

    [~,errors_ax_P1_noML] = convergenzaErrorePriori(...
        'P1',false,false,...
        true_sol_handle,grad_true_sol_handle,...
        'QuadratoMisto'...
        );

    noML_res(end+1,:,:) = errors_ax_P1_noML;

    [h_ax,errors_ax_P1_yesML] = convergenzaErrorePriori(...
        'P1',false,true,...
        true_sol_handle,grad_true_sol_handle,...
        'QuadratoMisto'...
        );

    yesML_res(end+1,:,:) = errors_ax_P1_yesML;
    toc
end

polyfit(log(h_ax), log(yesML_res(end,:,2)),1)

%%
figure(1)
subplot(1,2,1)
loglog(h_ax,noML_res(:,:,2)) 
title('senza ML')
xlabel('h')
ylabel('Errore L2')

subplot(1,2,2)
loglog(h_ax,yesML_res(:,:,2))
title('con ML')
xlabel('h')
ylabel('Errore L2')


sgtitle('Errore L2 FEM P1: ML vs no ML')

figure(2)
subplot(1,2,1)
loglog(h_ax,noML_res(:,:,3)) 
title('senza ML')
xlabel('h')
ylabel('Errore H1')

subplot(1,2,2)
loglog(h_ax,yesML_res(:,:,3))
title('con ML')
xlabel('h')
ylabel('Errore H1')


sgtitle('Errore H1 FEM P1: ML vs no ML')