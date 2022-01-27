clear all
close all
clc
if(~exist('assemblaEllittico'))
    addpath('Funzioni')
end

global geom
global problem

% Parametri per i P1
% Pk = 'P1';
% n_steps = 10;
% T = 4;

% Parametri per i P2
Pk = 'P2';
n_steps = 10;
T = 4;

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

area_ax = logspace(log10(0.1), log10(0.0001));
%area_ax = area_ax(1:5:end);

errori = [];

for area = area_ax
    find(area == area_ax)
    filename = append('Triangolazioni/QuadratoMisto/',num2str(area),'.mat');
    load(filename,'geom');
    if isequal(Pk,'P2') && (geom.nelements.nVertexes == length(geom.elements.coordinates))
        prepP2
    end
    
    [u,uD] = assemblaParabolico(Pk,T,n_steps,false,false,false);
    
    % Assembla la soluzione
    utilde = zeros(length(geom.pivot.pivot),n_steps+1);
    
    utilde(geom.pivot.pivot>0,:) = u;
    utilde(geom.pivot.pivot<0,:) = uD;
    
    errori(end+1) = calcolaErrorePrioriParabolico(true_sol_handle,utilde,Pk,T/n_steps);
end

polyfit(log(sqrt(area_ax)),log(errori),1)

%% Export to LATEX
%writematrix([sqrt(area_ax);errori]',"parabolico_convergenza_spazio_"+Pk+".csv")

%%
plot(log(sqrt(area_ax)),log(errori),'-o')
xlabel("log(h)")
ylabel("errore L2")
title("Convergenza errore L2 al variare di h per delta t = "+num2str(T/n_steps)+" con "+Pk)