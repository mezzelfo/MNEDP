clear all
close all
clc
if(~exist('assemblaEllittico'))
    addpath('Funzioni')
end


global geom;
global manifold;
global problem;


manifold.name = 'torus';
manifold.parameters.r1 = 1;
manifold.parameters.r2 = 2;

manifold_geodiff();
manifold_triangulator(100);
%%
problem.rho = @(x) 0*x(1,:)+1;
problem.epsilon = @(x) 0*x(1,:)+1;
problem.beta = @(x) 0*x;
problem.sigma = @(x) 0*x(1,:);
problem.f = @(x,t) 0*x(1,:);
%problem.iniziale = @(x) exp(-8*vecnorm(x-[pi;pi]).^2);
problem.iniziale = @(x) exp(-1./(1/2-vecnorm(x-[pi;pi]).^2)) .* (vecnorm(x-[pi;pi]).^2 < 1/2);
problem.v_iniziale = @(x) 0*x(1,:);

X = geom.elements.coordinates(:,1);
Y = geom.elements.coordinates(:,2);

T = 8;
n_passi = 200;
beta_Newmark = 1/4;
u = manifold_assemblaIperbolico('P1',T,n_passi,beta_Newmark);

utilde = u(geom.pivot.pivot,:);

%% Plot soluzione
map_surf = manifold.numeric.P(geom.elements.coordinates');
figure(1)
%pause
for t = 1:n_passi+1
    subplot(1,2,1)
    trisurf(geom.elements.triangles, map_surf(1,:), map_surf(2,:), map_surf(3,:), utilde(:,t))
    shading interp
    axis equal
    view(50,35)
    
    subplot(1,2,2)
    trisurf(geom.elements.triangles, geom.elements.coordinates(:,1), geom.elements.coordinates(:,2), utilde(:,t))
    view(0,90)
    axis([0 2*pi 0 2*pi])
    shading interp
    axis square
    
    pause(0.05)
end
return
%% Esporta come eps

tiledlayout(3,5,'TileSpacing','Compact','Padding','Compact')
for i = 1:15
    nexttile
    trisurf(geom.elements.triangles, map_surf(1,:), map_surf(2,:), map_surf(3,:), utilde(:,10*i))  
    shading interp
    axis equal
    %view(50,35)
    view(0,50)
    set(gca,'visible','off')
end


%% Export to Blender
% v_normale = manifold.numeric.normal(geom.elements.coordinates');
% writematrix(geom.elements.coordinates,'geom.elements.coordinates.csv');
% writematrix(geom.elements.triangles,'geom.elements.triangles.csv');
% writematrix(utilde,'utilde.csv');
% writematrix(v_normale','normal.csv');

% import bpy
% import numpy as np
% 
% uv_coords = np.genfromtxt('/home/luca/geom.elements.coordinates.csv',delimiter=',')
% r1 = 1
% r2 = 2
% coords = np.stack((
%     np.cos(uv_coords[:,0])*(r1*np.cos(uv_coords[:,1])+r2),
%     np.sin(uv_coords[:,0])*(r1*np.cos(uv_coords[:,1])+r2),
%     r1*np.sin(uv_coords[:,1])
% )).T
% 
% triangles = np.genfromtxt('/home/luca/geom.elements.triangles.csv',delimiter=',').astype(int)-1
% 
% animation = np.genfromtxt('/home/luca/utilde.csv',delimiter=',')
% 
% normals = np.genfromtxt('/home/luca/normal.csv',delimiter=',')
% 
% mesh = bpy.data.meshes.new("MATLABMesh")
% obj = bpy.data.objects.new(mesh.name, mesh)
% col = bpy.data.collections.get("Collection")
% col.objects.link(obj)
% mesh.from_pydata(coords.tolist(),[],triangles.tolist())
% 
% #ANIMATION
% sk_basis = obj.shape_key_add(name='Basis')
% sk_basis.interpolation = 'KEY_LINEAR'
% obj.data.shape_keys.use_relative = False
% 
% for frame in range(animation.shape[1]):
%     bpy.context.scene.frame_set(frame+1)
%     sk = obj.shape_key_add(name='Deform',from_mix=False)
%     sk.interpolation = 'KEY_LINEAR'
%     for i in range(len(coords)):
%         sk.data[i].co.x += 5*animation[i][frame]*normals[i][0]
%         sk.data[i].co.y += 5*animation[i][frame]*normals[i][1]
%         sk.data[i].co.z += 5*animation[i][frame]*normals[i][2]

