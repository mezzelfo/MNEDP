clear all
close all
clc
if(~exist('assemblaEllittico'))
    addpath('Funzioni')
end

%%
global geom
generaTriangolazione(0.0003, [0 0;1 0;1 1;0 1], [1 1 1 1], [3 3 3 3]);

%%
global problem


problem.beta = @(x) x*0;
problem.sigma = @(x) x(1,:)*0+1;
problem.f = @(x) x(1,:)*0+1;
problem.bordo_dirichlet = @(x, marker) x(1,:)*0;
problem.bordo_neumann = @(x, marker) x(1,:)*NaN;

figure(2)
clf
tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'none'); 

epsilon_ax = [1e-3, 1e-5, 1e-7];
for eps_idx = 1:3
    problem.epsilon = @(x) x(1,:)*0+epsilon_ax(eps_idx);
    
    nexttile(eps_idx)
    [A,b,A_dirichlet,u_dirichlet,b_neumann,~,~] = assemblaEllittico('P1',false,false,false);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    utilde = zeros(geom.nelements.nVertexes,1);
    utilde(geom.pivot.pivot>0) = u;
    utilde(geom.pivot.pivot<0) = u_dirichlet;
    trisurf(geom.elements.triangles,geom.elements.coordinates(:,1),geom.elements.coordinates(:,2),utilde)
    title("P1, no Mass Lumping, "+"epsilon="+num2str(epsilon_ax(eps_idx)))
    
    
    nexttile(eps_idx+3)
    [A,b,A_dirichlet,u_dirichlet,b_neumann,~,~] = assemblaEllittico('P1',false,true,false);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    utilde = zeros(geom.nelements.nVertexes,1);
    utilde(geom.pivot.pivot>0) = u;
    utilde(geom.pivot.pivot<0) = u_dirichlet;
    trisurf(geom.elements.triangles,geom.elements.coordinates(:,1),geom.elements.coordinates(:,2),utilde)
    title("P1, con Mass Lumping, "+"epsilon="+num2str(epsilon_ax(eps_idx)))
end


% prepP2();
% for eps_idx = 1:3
%     problem.epsilon = @(x) x(1,:)*0+epsilon_ax(eps_idx);
%     nexttile(3*eps_idx);
%     [A,b,A_dirichlet,u_dirichlet,b_neumann,~,~] = assemblaEllittico('P2',false,false,false);
%     u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
%     utilde = zeros(geom.nelements.nVertexes,1);
%     for v = 1:length(geom.pivot.pivot)
%         vv = geom.pivot.pivot(v);
%         if vv>0
%             utilde(v) = u(vv);
%         else
%             utilde(v) = u_dirichlet(-vv);
%         end
%     end
%     
%     [X,Y,Z] = evaluate_solution_at(utilde,[0 1 0 1],'P2');
%     surf(X,Y,Z)
% end


