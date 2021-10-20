function [Center,r2]=find_bordercircle_info (V1,V2);
%
% [xC,yC,r2] = find_circuminfo_of_memorized_tria (T)
%
% This function considers the bordercircle (*see attachments) of a segment
% determined by the coordinates of his verteces. The function calculates
% the center and the radius of the bordercircle
%
% Output & Input of this function are:
%
% [Center(1),Center(2)]: the coordinates of the bordercircle center
% r2: the square of the radius of the bordercircle
%
% [V1(1),V1(2)]: he coordinates of one of the defining vertex
% [V2(1),V2(2)]: he coordinates of the other defining vertex
%

global TV V

Center(1) = ( V1(1)+V2(1) )/2;
Center(2) = ( V1(2)+V2(2) )/2;
r2= (V1(1)-Center(1))^2 + (V1(2)-Center(2))^2;

return