close all
clear all
clc

Sample_Square_Dirichlet()





epsilon = 1;
betax = 0;
betay = 0;
gamma = 0;
f = @(P) 32*(P(1)*(1-P(1))+P(2)*(1-P(2)));

[A,b] = assembla(geom,epsilon,betax,betay,gamma,f);

u = A\b;

sol = zeros(length(geom.elements.coordinates),1);
sol(geom.pivot.pivot > 0) = u;

X = geom.elements.coordinates(:,1);
Y = geom.elements.coordinates(:,2);

subplot(2,2,1)
trisurf(geom.elements.triangles,X,Y,sol);

subplot(2,2,2)
trisurf(geom.elements.triangles,X,Y,16*X.*(1-X).*Y.*(1-Y));

subplot(2,2,3)
trisurf(geom.elements.triangles,X,Y,sol-16*X.*(1-X).*Y.*(1-Y));

