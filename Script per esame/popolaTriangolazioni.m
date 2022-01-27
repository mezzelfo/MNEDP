clear all
close all
clc
if(~exist('assemblaEllittico'))
    addpath('Funzioni')
end

%QuadratoMisto
InputVertex = [0 0; 1 0; 1 1; 0 1];
InputVertexValues = [1 1 1 1];
BoundaryValues = [3 4 5 7];

% %NonConvesso
% InputVertex = [0 0; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1 0];
% BoundaryValues = [1 1 1 1 1 1 1 1];
% InputVertexValues = [1 1 1 1 1 1 1 1];



area_ax = logspace(log10(0.1), log10(0.0001));
for i = 1:length(area_ax)
    generaTriangolazione(area_ax(i),InputVertex,InputVertexValues,BoundaryValues,'QuadratoMisto')
end