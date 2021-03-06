close all
clear all
clc

Sample_Square_Dirichlet()

epsilon = @(x) 1;
betax = @(x) 0;
betay = @(x) 0;
gamma = @(x) 0;
f = @(P) 32*(P(1)*(1-P(1))+P(2)*(1-P(2)));
bordo_dirichlet = @(P, type) 5;

[A,b,A_dirichlet,u_dirichlet] = assembla(geom,epsilon,betax,betay,gamma,f, bordo_dirichlet);

u = A\(b-A_dirichlet*u_dirichlet);

sol = zeros(length(geom.elements.coordinates),1);
sol(geom.pivot.pivot > 0) = u;
sol(geom.pivot.pivot < 0) = u_dirichlet;

X = geom.elements.coordinates(:,1);
Y = geom.elements.coordinates(:,2);

subplot(2,2,1)
trisurf(geom.elements.triangles,X,Y,sol);

subplot(2,2,2)
trisurf(geom.elements.triangles,X,Y,5+16*X.*(1-X).*Y.*(1-Y));

subplot(2,2,3)
trisurf(geom.elements.triangles,X,Y,sol-5-16*X.*(1-X).*Y.*(1-Y));

