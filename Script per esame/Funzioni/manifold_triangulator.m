function manifold_triangulator(N)
global geom;
global manifold;

% Assumiamo che il dominio sia il rettangolo descritto in manifold.uv_bounds
X=linspace(manifold.uv_bounds(1,1),manifold.uv_bounds(1,2),N)';
Y=linspace(manifold.uv_bounds(2,1),manifold.uv_bounds(2,2),N)';
[X,Y] = meshgrid(X,Y);
X = X(:);
Y = Y(:);
TR = delaunayTriangulation(X,Y);

% Assumiamo:
% 1) x punto interno rettangolo => P(x) punto interno varietà
% 2) x punto di bordo rettangolo t.c. non esiste y t.c. P(x) = P(y) =>
% necessarie condizione di Dirichlet per x
% 3) x punto di bordo rettangolo t.c. esiste y t.c. P(x) = P(y) =>
% P(x) punto interno di varietà

mapped = manifold.numeric.P(TR.Points');

mapped_distance = squareform(pdist(mapped'));
mapped_near = mapped_distance < 100*eps;
[~,rappresentante_classe] = max(mapped_near);

% rappresentante_classe = @(e) e;
tol = 100*eps;

% +---3---+
% |       |
% 4       2
% |       |
% +---1---+

lato3 = abs(Y-manifold.uv_bounds(2,2)) < tol;
lato1 = abs(Y-manifold.uv_bounds(2,1)) < tol;
lato2 = abs(X-manifold.uv_bounds(1,2)) < tol;
lato4 = abs(X-manifold.uv_bounds(1,1)) < tol;

lati = [lato1 lato2 lato3 lato4];

Ne = length(mapped);
pivot = -ones(Ne,1);
i = 1;
z = 1;
for e = 1:Ne
    if any(lati(e,:) .* manifold.dirichlet_borders)
        %dirichlet
        incollato = rappresentante_classe(e);
        if incollato == e
            pivot(e) = -z;
            z = z+1;
        else
            pivot(e) = pivot(incollato);
        end
    else
        incollato = rappresentante_classe(e);
        if incollato == e
            pivot(e) = i;
            i = i+1;
        else
            pivot(e) = pivot(incollato);
        end
    end
%     if bordo_rettangolo(e) == 0 %Tutti i punti interni devono essere liberi
%         pivot(e) = i;
%         i = i + 1;
%     elseif cardinalita_classe(e) == 1 %Punti non incollati sul bordo del rettangolo => bordo della varietà
%         pivot(e) = -z;
%         z = z+1;
%     else %Punti di bordo del rettangolo che devono essere incollati
%         incollato = rappresentante_classe(e);
%         if incollato == e
%             pivot(e) = i;
%             i = i+1;
%         else
%             pivot(e) = pivot(incollato);
%         end
%     end
end

for t = 1:length(TR.ConnectivityList)
    idx = TR.ConnectivityList(t,:);
    coords = TR.Points(idx,:);
    a = polyarea(coords(:,1),coords(:,2));
    geom.support.TInfo(t).Area = a;
end

geom.pivot.pivot = pivot;
geom.elements.triangles = TR.ConnectivityList;
geom.elements.coordinates = TR.Points;
geom.nelements.nTriangles = size(TR.ConnectivityList);
geom.nelements.nVertexes = length(pivot);
end