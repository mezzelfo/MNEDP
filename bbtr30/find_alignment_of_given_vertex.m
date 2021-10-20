function [Position] = find_alignment_of_given_vertex(V1,V2,V3)
%
% [Position] =  = find_alignment_of_given_vertex(V1,V2,V3)
%
% This function checks if verteces V1, V2 & V3 are aligned or not. The
% function calculate a quantity related to the area of triangle V1-V2-V3
% using relation
%
% Area = 1/2 det [x1 x2 x3; y1 y2 y3; 1 1 1]
%
% The function uses global tolerance
%
% Output & Input of this function are:
% Position: a flag that assumes value 0 if verteces are aligned and 1 if
%           not
% V1,V2 : two verteces to control given by reference
% V[3(1) V3(2)]: a vertex given by coordinates


global V
global tolerance

Position = ( V(V2,2)*( V(V1,1)-V3(1) ) + V3(2)*( V(V2,1)-V(V1,1) )  + V(V1,2)*( V3(1)-V(V2,1) ) );

if ( (Position < tolerance) & (Position > -tolerance) );
    Position = 0;
else
    Position = 1;
end

return



    

