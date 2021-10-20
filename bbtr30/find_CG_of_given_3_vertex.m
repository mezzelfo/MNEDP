function [CG]=find_CG_of_given_3_vertex (V1,V2,V3);
%
% [CG]=find_CG_of_given_3_vertex (V1,V2,V3);
%
% This function find the center of gravity of three verteces, calculating
% the arithmetic average
%
% Output & Input of this function are:
%
% V1: a 2*1 array containing the coordinates of first vertex;
% V2: a 2*1 array containing the coordinates of second vertex;
% V3: a 2*1 array containing the coordinates of third vertex;
%
% CG: a 2*1 array containing the coordinates of center of gravity;
%

CG = [ (V1(1)+V2(1)+V3(1))/3 (V1(2)+V2(2)+V3(2))/3 ];

return