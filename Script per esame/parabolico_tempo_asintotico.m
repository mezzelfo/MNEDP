clear all
%close all
clc
if(~exist('assemblaEllittico'))
    addpath('Funzioni')
end

global geom
global problem

% Parametri per i P1
Pk = 'P1';
%area = 0.001;
area = 0.0001;
T = 4;

% Parametri per i P2
% Pk = 'P2';
% area = 0.001;
% T = 4;

% Genera una triangolazione con h piccolo
generaTriangolazione(area, [0 0;1 0;1 1;0 1], [1 1 1 1], [1 2 1 1])

true_sol_handle = @(x,t) 100*exp(-t).*(x(1,:).^2+x(2,:).^2);

problem.epsilon = @(x) x(1,:)*0+1;
problem.beta = @(x) 0*x;
problem.sigma = @(x) 0*x(1,:);
problem.f = @(x,t) -100*exp(-t).*(4+x(1,:).^2+x(2,:).^2);
problem.bordo_dirichlet = @(x,t, marker) true_sol_handle(x,t);
problem.bordo_neumann = @(x,t, marker) 100*2*exp(-t).*x(1,:);
problem.rho = @(x) 0*x(1,:) + 1;
problem.iniziale = @(x) true_sol_handle(x,0);

%% Calcola errori
%n_steps_ax = floor(logspace(0.4,3));
n_steps_ax = round(logspace(1,3));
%n_steps_ax = [n_steps_ax(1:5),n_steps_ax(6:5:end)];


errori = [];%zeros(1,length(n_steps_ax));

if isequal(Pk,'P2') && (geom.nelements.nVertexes == length(geom.elements.coordinates))
    prepP2
end

tic
for n_steps = n_steps_ax
    n_steps
    [u,uD] = assemblaParabolico(Pk,T,n_steps,false,false,false);
    
    % Assembla la soluzione
    utilde = zeros(length(geom.pivot.pivot),n_steps+1);
    
    utilde(geom.pivot.pivot>0,:) = u;
    utilde(geom.pivot.pivot<0,:) = uD;
    %     for t = 1:n_steps+1
    %         for v = 1:length(geom.pivot.pivot)
    %             ii = geom.pivot.pivot(v);
    %             if ii > 0 % Dof
    %                 utilde(v,t) = u(ii,t);
    %             else
    %                 utilde(v,t) = uD(-ii,t);
    %             end
    %         end
    %     end
    errori(end+1) = calcolaErrorePrioriParabolico(true_sol_handle,utilde,Pk,T/n_steps);
end
toc

polyfit(log(T./n_steps_ax),log(errori),1)

%%
loglog(T./n_steps_ax,errori,'o')
%%
%writematrix([T./n_steps_ax;errori]',"parabolico_convergenza_tempo_"+Pk+".csv")

%%
% trisurf(geom.elements.triangles(:,1:3),...
% geom.elements.coordinates(:,1),...
% geom.elements.coordinates(:,2),...
% utilde(:,end))

%%
% figure(1)
% plot(log(T./n_steps_ax),log(errori),'o-')
% xlabel("log(delta t)")
% ylabel("errore L2")
% title("Convergenza errore L2 al variare di delta t per h = "+num2str(sqrt(area))+" con "+Pk)
% 
% %%
% clear all
% clc
%p1 = readmatrix('parabolico_convergenza_tempo_.csv');
%p2 = readmatrix('parabolico_convergenza_tempo_P2.csv');
