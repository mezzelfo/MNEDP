function assembleManifold()
% Assuming a rectangular chart with borders labeld according to the following scheme
% +---3---+
% |       |
% 4       2
% |       |
% +---1---+
% the function takes a global struct "manifold" containing
% - uv_bounds: bounds of the rectangular chart
% - dirichlet_borders: a binary array with 1 in position i iff border i
%       is a boundary with Dirichlet conditions
% - P: the symbolic parametrization map possibly dependent on some
%       parameters
% and computes the following useful quantities
% - symbolic.g: symbolic metric tensor
% - symbolic.ginv: symbolic inverse metric tensor
% - symbolic.sqrtdetg: symbolic volume form
% - symbolic.normal: symbolic normal field
% - numeric.P: parametrization's function handle
% - numeric.ginv: inverse metric tensor's function handle
% - numeric.sqrtdetg: volume form's function handle
% - numeric.normal: normal field's function handle
% of a 2-manifold embedded in R^3 with standard metric.
% If the function is called with a manifold.name value in the following list
% {torus, klein_bottle, sphere, cylinder_x, cylinder_y, spherical_cap} then
% the correct uv_bounds, dirichlet_borders and P are already coded.

global manifold;
syms P(u,v)
assume([u,v],'real')


% Pre-coded surfaces
switch manifold.name
    case 'torus'
        r1 = manifold.parameters.r1;
        r2 = manifold.parameters.r2;
        P(u,v) = [...
            cos(u)*(r1*cos(v)+r2);...
            sin(u)*(r1*cos(v)+r2);...
            r1*sin(v)];
        manifold.uv_bounds = [0,2*pi;0,2*pi];
        manifold.dirichlet_borders = [0 0 0 0];
    case 'klein_bottle'
        r = manifold.parameters.r;
        P(u,v) = [...
            (r+cos(u/2)*sin(v)-sin(u/2)*sin(2*v))*cos(u);...
            (r+cos(u/2)*sin(v)-sin(u/2)*sin(2*v))*sin(u);...
            sin(u/2)*sin(v) + cos(u/2)*sin(2*v)];
        manifold.uv_bounds = [0,2*pi;0,2*pi];
        manifold.dirichlet_borders = [0 0 0 0];
    case 'sphere'
        r = manifold.parameters.r;
        P(u,v) = [...
            r*cos(u)*sin(v);...
            r*sin(u)*sin(v);...
            r*cos(v)];
        manifold.uv_bounds = [0,2*pi;0,pi];
        manifold.dirichlet_borders = [0 0 0 0];
    case 'cylinder_x'
        r = manifold.parameters.r;
        P(u,v) = [...
            r*cos(u);...
            r*sin(u);...
            v];
        manifold.uv_bounds = [0,2*pi;0,2*pi];
        manifold.dirichlet_borders = [1 0 1 0];
    case 'cylinder_y'
        r = manifold.parameters.r;
        P(u,v) = [...
            r*cos(v);...
            r*sin(v);...
            u];
        manifold.uv_bounds = [0,2*pi;0,2*pi];
        manifold.dirichlet_borders = [0 1 0 1];
    case 'spherical_cap'
        r = manifold.parameters.r;
        P(u,v) = [...
            u;...
            v;...
            r^2-u^2-v^2];
        manifold.uv_bounds = [-r/2,r/2;-r/2,r/2];
        manifold.dirichlet_borders = [1 1 1 1];
    case 'paraboloid'
        r = manifold.parameters.r;
        P(u,v) = [...
            u;...
            v;...
            u^2+v^2];
        manifold.uv_bounds = [-r/2,r/2;-r/2,r/2];
        manifold.dirichlet_borders = [1 1 1 1];
    case 'mobius_strip'
        r = manifold.parameters.r;
        w = manifold.parameters.w;
        P(u,v) = [...
            (r+u*cos(v/2))*cos(v);...
            (r+u*cos(v/2))*sin(v);...
            u*sin(v/2)];
        manifold.uv_bounds = [-w,w;0,2*pi];
        manifold.dirichlet_borders = [0,1,0,1];
    otherwise
        disp('No manifold name specified.')
end

% Tangent vector fields
Pu = diff(P,u);
Pv = diff(P,v);

J = [Pu, Pv];
g = simplify(J'*J); % Metric tensor
ginv = simplify(inv(g)); % Inverse metric tensor
sqrtdetg = simplify(sqrt(abs(det(g)))); % Volume form
N = simplify(cross(Pu,Pv)/norm(cross(Pu,Pv))); % Normal map

manifold.symbolic.P = P;
manifold.symbolic.g = g;
manifold.symbolic.ginv = ginv;
manifold.symbolic.sqrtdetg = sqrtdetg;
manifold.symbolic.normal = N;

manifold.numeric.P = matlabFunction(P,'Vars',{[u;v]});
manifold.numeric.ginv = matlabFunction(ginv,'Vars',{[u;v]});
manifold.numeric.sqrtdetg = matlabFunction(sqrtdetg,'Vars',{[u;v]});
manifold.numeric.normal = matlabFunction(N,'Vars',{[u;v]});
end
