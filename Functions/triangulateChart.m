function triangulateChart(N,M)
% Given a pre-assembled "manifold" struct this produces a struct "geom"
% containing a triangulation based on a NxM grid of the domain specified by
% manifold.uv_bounds. The degrees of freedom are computed to glue together
% points of the chart which are identified by the parametrization map
% manifold.numeric.P.
% geom contains
% pivot.pivot: degrees of freedom indices, negative if dirichlet points
% elements.triangles: list of the triangles and their vertices' indices
% elements.coordinates: vertices uv coordinates
% nelements.nTriangles: number of triangles
% nelements.nVertexes: number of vertices
%geom.support.TInfo(t).Area: area of triangle t

global geom;
global manifold;

X=linspace(manifold.uv_bounds(1,1),manifold.uv_bounds(1,2),N)';
Y=linspace(manifold.uv_bounds(2,1),manifold.uv_bounds(2,2),M)';
[X,Y] = meshgrid(X,Y);
X = X(:);
Y = Y(:);
TR = delaunayTriangulation(X,Y); % Triangulation

% Assume
% 1) x interior point of the rectangle => P(x) interior point of the
%   manifold
% 2) x boundary point of the chart such that there is no y with P(x) = P(y) =>
%   Dirichlet conditions for x are needed (x is not glued to any other point)
% 3) x boundary point of the chart such that exists y with P(x) = P(y) =>
%   P(x) interior point of the manifold (x and y are glued)

mapped = manifold.numeric.P(TR.Points'); % Chart points mapped in R^3

mapped_distance = squareform(pdist(mapped')); % Distance matrix
mapped_near = mapped_distance < 100*eps;
% Index of the representative point for each class of glued points
[~,class_representative] = max(mapped_near);

% rappresentante_classe = @(e) e;
tol = 100*eps;

% +---3---+
% |       |
% 4       2
% |       |
% +---1---+

% Side points indices masks
side1 = abs(Y-manifold.uv_bounds(2,1)) < tol;
side2 = abs(X-manifold.uv_bounds(1,2)) < tol;
side3 = abs(Y-manifold.uv_bounds(2,2)) < tol;
side4 = abs(X-manifold.uv_bounds(1,1)) < tol;

sides = [side1 side2 side3 side4];

Ne = length(mapped); % Number of points
pivot = -ones(Ne,1); % degrees of freedom array
i = 1;
z = 1;
for e = 1:Ne
    if any(sides(e,:) .* manifold.dirichlet_borders)
        % The point is on a Dirchlet border
        glued = class_representative(e);
        if glued == e
            pivot(e) = -z;
            z = z+1;
        else
            pivot(e) = pivot(glued);
        end
    else
        % The point is in the manifold's interior
        glued = class_representative(e);
        if glued == e
            pivot(e) = i;
            i = i+1;
        else
            pivot(e) = pivot(glued);
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

% Precompute the triangles' areas
for t = 1:length(TR.ConnectivityList)
    idx = TR.ConnectivityList(t,:);
    coords = TR.Points(idx,:);
    a = polyarea(coords(:,1),coords(:,2));
    geom.support.TInfo(t).Area = a;
end


geom.pivot.pivot = pivot; % Degrees of freedom indices
geom.elements.triangles = TR.ConnectivityList; % Triangles list
geom.elements.coordinates = TR.Points; % Vertices coordinates list
geom.nelements.nTriangles = size(TR.ConnectivityList); % Number of triangles
geom.nelements.nVertexes = length(pivot); % Number of vertices
end