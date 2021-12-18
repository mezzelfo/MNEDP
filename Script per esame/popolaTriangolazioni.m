clear all
close all
clc
if(~exist('assemblaEllittico'))
     addpath('Funzioni')
end


%InputVertex = [0 0; 1 0; 0.3 0.5; 1 1; 0 1];
%BoundaryValues = 0*InputVertex+1;

InputVertex = [0 0; 0.57 0.15; 1 0; 0.74 0.5; 1 1; 0.54 0.83; 0 1; 0.18 0.52];
BoundaryValues = 0*InputVertex+1;


area_ax = logspace(log10(0.1), log10(0.0001));
%area_ax = area_ax(1:10:end);
for i = 1:length(area_ax)
    generaTriangolazione(area_ax(i),InputVertex,BoundaryValues,'Prova')
end