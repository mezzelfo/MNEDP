function [RemovedTriangle,TriangleKind] = check_if_triangle_is_to_force_in_grid (T,RemovedTriangle,TriangleKind);
%
% [RemovedTriangle,TriangleKind] = check_if_triangle_is_to_force_in_grid (T,RemovedTriangle,TriangleKind);
%
% This function considers a triangle T that array TriangleKind says if
% belongs or not to the domain. If T lyes outside the domain the 
% function forces T in grid if a small angle 
% (< than global parameter LimAngleToForceInGrid) is present in T and if it is
% bounded by triangles beloning or forced into domain. The function then checks the
% neighbour triangle opposite to the small angle, and eventually forces
% also if (the function calls itself)
%
% Output & Input of this function are:
%
% RemovedTriangles: a boolean that becomes true
%   if some exterior triangle is removed from triangulation
%
% TriangleKind: is an array builded in the way that 
%   TriangleKind(i) = 0 if i doesn't belong to domain and
%   TriangleKind(i) = 1 if i does
%
% T is the triangle to check
%

global LimAngleToForceInGrid
global TT

% Create Recognizer array, that says if TT(T,I) must be considered
for I = 1:3
    if TT(T,I) == -1
        Recognizer(I) = 0;
    elseif TriangleKind(TT(T,I)) == 0
        Recognizer(I) = 0;
    else
        Recognizer(I) = 1;
    end
end

if TriangleKind(T) == 0 
    
    % This is a triangle to delete. Search if it has two (domain border o
    if Recognizer(1) == 1

        if Recognizer(2) == 1
            
            Angle = evaluate_angle (T,TT(T,6),TT(T,4),TT(T,5));
            if Angle < LimAngleToForceInGrid
                TriangleKind(T) = 1;
                RemovedTriangle = true;
                [RemovedTriangle TriangleKind] = check_if_triangle_is_to_force_in_grid (TT(T,3),RemovedTriangle,TriangleKind);
            elseif Recognizer(3) == 1
 
                % Unfortunatly the triangle is fully bounded by domain
                % limits. Three checks must be made.
                Angle = evaluate_angle (T,TT(T,5),TT(T,4),TT(T,6));
                if Angle < LimAngleToForceInGrid
                    TriangleKind(T) = 1;
                    RemovedTriangle = true;
                    [RemovedTriangle TriangleKind] = check_if_triangle_is_to_force_in_grid (TT(T,2),RemovedTriangle,TriangleKind);
                else
                    Angle = evaluate_angle (T,TT(T,4),TT(T,5),TT(T,6));
                    if Angle < LimAngleToForceInGrid
                        TriangleKind(T) = 1;
                        RemovedTriangle = true;
                        [RemovedTriangle TriangleKind] = check_if_triangle_is_to_force_in_grid (TT(T,3),RemovedTriangle,TriangleKind);
                    end
                end
                
            end
            
        elseif Recognizer(3) == 1
            
            Angle = evaluate_angle (T,TT(T,5),TT(T,4),TT(T,6));
            if Angle < LimAngleToForceInGrid
                TriangleKind(T) = 1;
                RemovedTriangle = true;
                [RemovedTriangle TriangleKind] = check_if_triangle_is_to_force_in_grid (TT(T,2),RemovedTriangle,TriangleKind);
            end
            
        end
        
    elseif (Recognizer(2) == 1) & (Recognizer(3) == 1)
        
        Angle = evaluate_angle (T,TT(T,4),TT(T,5),TT(T,6));
        if Angle < LimAngleToForceInGrid
            TriangleKind(T) = 1;
            RemovedTriangle = true;
            [RemovedTriangle TriangleKind] = check_if_triangle_is_to_force_in_grid (TT(T,3),RemovedTriangle,TriangleKind);
        end
        
    end
    
end

return