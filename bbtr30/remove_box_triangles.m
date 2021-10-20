function [RemovedTriangles] = remove_box_triangles();
%
% [RemovedTriangles] = function remove_box_triangles();
%
% This function considers the triangulation and, using Autoidentify
% Algorithm, determines which triangles must be removed. Then the function
% checks if some triangle must be included in the grid (RemoveTriangles in
% this case become true) even if outside from domain (*see attachments),
% and removes from grid the triangle that must be deletes
%
% Output of this function is
%
% RemovedTriangles: a boolean that becomes true if some exterior triangle
%       is removed from triangulation
%

global nT

% Identify different regions
TriangleKind = autoidentify_grid; %0 delete, 1 keep

% Set RemovedTriangles
RemovedTriangles = false;

% Determine if some triangle should not be deleted due to it's a small
% angle
for iT = 1 : nT
    [RemovedTriangles TriangleKind] = check_if_triangle_is_to_force_in_grid (iT,RemovedTriangles,TriangleKind);    
end

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