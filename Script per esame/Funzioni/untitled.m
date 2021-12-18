InputVertex = [0 0; 1 0; 1 1; 0 1];
BoundaryValues = [3 4 5 7];

area_ax = logspace(log10(0.1), log10(0.0001));
%area_ax = area_ax(1:10:end);
h_ax = sqrt(area_ax);
errors_ax = zeros(length(area_ax),3);
for i = 1:length(area_ax)
    generaTriangolazione(area_ax(i),InputVertex,BoundaryValues,'QuadratoMisto')
end