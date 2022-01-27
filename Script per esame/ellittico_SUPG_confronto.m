clear all
close all
clc
if(~exist('assemblaEllittico'))
    addpath('Funzioni')
end

%%
global geom
%generaTriangolazione(0.0003, [0 0;1 0;1 1;0 1], [1 1 1 1], [3 3 3 3])
generaTriangolazione(0.001, [0 0;1 0;1 1;0 1], [1 1 1 1], [3 3 3 3])

%% %TODO: Chiedere a Berrone se bisogna tenere u_true fissa o f fissa.
global problem


problem.beta = @(x) x*0+1;
problem.sigma = @(x) x(1,:)*0;
problem.f = @(x) x(1,:)*0+1;
problem.bordo_dirichlet = @(x, marker) x(1,:)*0;
problem.bordo_neumann = @(x, marker) x(1,:)*NaN;


figure(1)
clf
tiledlayout(3,4, 'Padding', 'none', 'TileSpacing', 'none');

epsilon_ax = [1e-3, 1e-5, 1e-7];
for eps_idx = 1:3
    problem.epsilon = @(x) x(1,:)*0+epsilon_ax(eps_idx);
    
    nexttile(4*eps_idx-3);
    [A,b,A_dirichlet,u_dirichlet,b_neumann,~,~] = assemblaEllittico('P1',false,false,false);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    utilde = zeros(geom.nelements.nVertexes,1);
    utilde(geom.pivot.pivot>0) = u;
    utilde(geom.pivot.pivot<0) = u_dirichlet;
    trisurf(geom.elements.triangles,geom.elements.coordinates(:,1),geom.elements.coordinates(:,2),utilde)
    view(60,50)
    title("P1, no SUPG, "+"epsilon="+num2str(epsilon_ax(eps_idx)))
    
    
    nexttile(4*eps_idx-2);
    [A,b,A_dirichlet,u_dirichlet,b_neumann,~,~] = assemblaEllittico('P1',true,false,false);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    utilde = zeros(geom.nelements.nVertexes,1);
    utilde(geom.pivot.pivot>0) = u;
    utilde(geom.pivot.pivot<0) = u_dirichlet;
    trisurf(geom.elements.triangles,geom.elements.coordinates(:,1),geom.elements.coordinates(:,2),utilde)
    view(60,50)
    title("P1, con SUPG, "+"epsilon="+num2str(epsilon_ax(eps_idx)))
end

prepP2();
for eps_idx = 1:3
    problem.epsilon = @(x) x(1,:)*0+epsilon_ax(eps_idx);
    nexttile(4*eps_idx-1);
    [A,b,A_dirichlet,u_dirichlet,b_neumann,~,~] = assemblaEllittico('P2',false,false,false);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    utilde(geom.pivot.pivot>0) = u;
    utilde(geom.pivot.pivot<0) = u_dirichlet;
    [X,Y,Z] = evaluate_solution_at(utilde,[0 1 0 1],'P2');
    surf(X,Y,Z);
    view(60,50)
    title("P2, no SUPG "+"epsilon="+num2str(epsilon_ax(eps_idx)))
    
    nexttile(4*eps_idx);
    [A,b,A_dirichlet,u_dirichlet,b_neumann,~,~] = assemblaEllittico('P2',true,false,false);
    u = A\(b-A_dirichlet*u_dirichlet+b_neumann);
    utilde(geom.pivot.pivot>0) = u;
    utilde(geom.pivot.pivot<0) = u_dirichlet;
    [X,Y,Z] = evaluate_solution_at(utilde,[0 1 0 1],'P2');
    surf(X,Y,Z);
    view(60,50)
    title("P2, con SUPG "+"epsilon="+num2str(epsilon_ax(eps_idx)))
end

%%
% figure(2)
% clf
% triplot(triangulation(geom.elements.triangles(:,1:3),geom.elements.coordinates))
% hold on
% scatter(X(:),Y(:),3,'*')
% hold off