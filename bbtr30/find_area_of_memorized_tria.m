function [Area] = find_area_of_memorized_tria (T)
%
% [Area] = find_area_of_memorized_tria(T)
%
% This function calculates the area of a triangle, given the reference
% inthe global variable of the considered triangle. The function also
% checks if verteces are in anticlockwise order and eventually invert their
% order. The used way is to calculate
%
%
% 1/2 det [x1 x2 x3; y1 y2 y3; 1 1 1]
%
% Area: the area to be calculated;
% T: the reference to the chosen triangle

global TV V TT

Area = 1/2*( V(TV(T,2),2)*( V(TV(T,1),1)-V(TV(T,3),1) ) + V(TV(T,3),2)*( V(TV(T,2),1)-V(TV(T,1),1) )  + V(TV(T,1),2)*( V(TV(T,3),1)-V(TV(T,2),1) ) );

if Area < 0
    Area = -Area;
    
    temp = TV(T,1);
    TV(T,1) = TV(T,2);
    TV(T,2) = temp;
    
end

return



    

