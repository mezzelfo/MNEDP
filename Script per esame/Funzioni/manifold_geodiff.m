function manifold_geodiff()
global manifold;
syms P(u,v)
assume([u,v],'real')

% +---3---+
% |       |
% 4       2
% |       |
% +---1---+

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
    otherwise
        error('Unexpected manifold type. No geometry created.')
end

Pu = diff(P,u);
Pv = diff(P,v);
J = [Pu, Pv];
g = simplify(J'*J);
ginv = simplify(inv(g));
sqrtdetg = simplify(sqrt(abs(det(g))));
N = simplify(cross(Pu,Pv)/norm(cross(Pu,Pv)));

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
