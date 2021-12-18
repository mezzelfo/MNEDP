clear all
close all
clc
if(~exist('assemblaEllittico'))
     addpath('Funzioni')
end

%%
global problem
% true_sol_handle = @(x) sin(3*x(1,:)).*cos(4*x(2,:));
% grad_true_sol_handle = @(x) [3*cos(3*x(1,:)).*cos(4*x(2,:)); -4*sin(3*x(1,:)).*sin(4*x(2,:))];
% 
% problem.epsilon = @(x) x(1,:)*0+1;
% problem.beta = @(x) [x(2,:);-x(1,:)];
% problem.sigma = @(x) x(1,:);
% problem.f = @(x) (25+x(1,:)).*true_sol_handle(x)+dot(problem.beta(x), grad_true_sol_handle(x));
% problem.bordo_dirichlet = @(x, marker) true_sol_handle(x);
% problem.bordo_neumann = @(x, marker) 3*cos(3*x(1,:)).*cos(4*x(2,:));

true_sol_handle = @(x) 5+16*x(1,:).*(1-x(1,:)).*x(2,:).*(1-x(2,:));
grad_true_sol_handle = @(x) [   16*x(2,:).*(1-x(2,:)).*(1-2*x(1,:));
    16*x(1,:).*(1-x(1,:)).*(1-2*x(2,:))];
problem.epsilon = @(x) x(1,:)*0 + 1;
problem.beta = @(x) x*0+[0;0];
problem.sigma = @(x) x(1,:)*0 + 1;
problem.bordo_dirichlet = @(x, marker) true_sol_handle(x);
problem.bordo_neumann = @(x,marker) 0;
problem.f = @(x) 32*(x(1,:).*(1-x(1,:))+x(2,:).*(1-x(2,:))) + problem.sigma(x).*true_sol_handle(x);
%% PLOT?
% model = fitlm(log(h_ax), log(errors_ax_P1(:,2)));
% h = model.plot;
% h(4).Color = 'green';
% h

%% P1 vs P2
[h_ax,errors_ax_P1] = convergenzaErrorePriori(...
    'P1',false,false,...
    true_sol_handle,grad_true_sol_handle,...
    'Prova'...
    );

% [h_ax,errors_ax_P2] = convergenzaErrorePriori(...
%     'P2',false,false,...
%     true_sol_handle,grad_true_sol_handle,...
%     'NonConvesso'...
%     );

% loglog(h_ax,errors_ax_P1(:,2),'-o',h_ax,errors_ax_P2(:,2),'-o')
% legend('P1','P2');

polyfit(log(h_ax), log(errors_ax_P1(:,1)),1)
polyfit(log(h_ax), log(errors_ax_P1(:,2)),1)
polyfit(log(h_ax), log(errors_ax_P1(:,3)),1)

% polyfit(log(h_ax), log(errors_ax_P2(:,1)),1)
% polyfit(log(h_ax), log(errors_ax_P2(:,2)),1)
% polyfit(log(h_ax), log(errors_ax_P2(:,3)),1)
%%
load("Triangolazioni\NonConvesso\0.00020236.mat");
[A,b,A_dirichlet,u_dirichlet,b_neumann,s,se] = assemblaEllittico('P1',false,false,false);
u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
utilde = zeros(geom.nelements.nVertexes,1);
utilde(geom.pivot.pivot>0) = u;
utilde(geom.pivot.pivot<0) = u_dirichlet;
trisurf(geom.elements.triangles,geom.elements.coordinates(:,1),geom.elements.coordinates(:,2),utilde)
