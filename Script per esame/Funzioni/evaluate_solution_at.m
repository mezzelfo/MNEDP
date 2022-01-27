function [X,Y,Z] = evaluate_solution_at(utilde,limits,Pk)
global geom
x_min = limits(1);
x_max = limits(2);
y_min = limits(3);
y_max = limits(4);
x = linspace(x_min,x_max,100);
y = linspace(y_min,y_max,100);
[X,Y] = meshgrid(x,y);
all_points = [X(:)'; Y(:)'];
Z = zeros(1,length(all_points));
Z_used = zeros(1,length(all_points));

switch Pk
    case 'P1'
        phi = {@(x) x(1,:); @(x) x(2,:); @(x) (1-x(1,:)-x(2,:))};
    case 'P2'
        phi = {@(x) 2.*x(1,:).*(x(1,:)-0.5); @(x) 2.*x(2,:).*(x(2,:)-0.5); @(x) 2.*(1-x(1,:)-x(2,:)).*((1-x(1,:)-x(2,:))-0.5);@(x) 4*x(1,:).*x(2,:);@(x) 4.*x(2,:).*(1-x(1,:)-x(2,:));@(x) 4.*x(1,:).*(1-x(1,:)-x(2,:))};
    otherwise
        error('Funzione implementata solo per P1 e P2');
end


for e=1:geom.nelements.nTriangles
    dof = geom.elements.triangles(e,:);
    coords = geom.elements.coordinates(dof,:)';
    d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
    %B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
    %pts = B \ (all_points-coords(:,3));
    area = geom.support.TInfo(e).Area;
    B_inv = [d(2,1) d(1,1); d(2,2) d(1,2)]/(2*area);
        
    pts = B_inv * (all_points-coords(:,3));
    selection = (pts(1,:) >= 0) & (pts(2,:) >= 0) & ((pts(1,:)+pts(2,:)) <= 1);
    selection = selection & not(Z_used);
    for j=1:length(dof)
        f = phi{j};
        Z = Z + utilde(dof(j))*f(pts).*selection;
    end
    Z_used = Z_used | selection;
end
Z = reshape(Z, length(x), length(y));
end

