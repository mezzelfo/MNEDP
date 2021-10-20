function [] = remove_outside_triangles();
%
% [] = function remove_outside_triangles();
%
% This function considers the triangulation and, using Autoidentify
% Algorithm, determines which triangles must be removed. 
%


global nT

% Identify different regions
TriangleKind = autoidentify_grid; %0 delete, 1 keep

% Create a plain list of triangles to be deleted
nTriangleToBeDeleted = 0;
for iT = 1 : nT
    
    if TriangleKind(iT) == 0
        nTriangleToBeDeleted = nTriangleToBeDeleted + 1;
        TriangleToBeDeleted(nTriangleToBeDeleted) = iT;
    end    
    
end


% Delete triangle that must be deleted
delete_triangles (TriangleToBeDeleted,nTriangleToBeDeleted);

return