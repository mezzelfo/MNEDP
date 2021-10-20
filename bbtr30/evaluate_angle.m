function [Angle] = evaluate_angle (T,OppositeB,AdjacentB1,AdjacentB2);
%
% Angle = evaluate_angle (T,OppositeB)
%
% This function calulates the angle of a memorized triangle T. The way used
% to recognize the angle among the three angles is to consider which are
% Borders in the triangle. To make the calculation we considered the
% relation:
%
%     1/2 A B sen gamma = AreaOfTriangle, 
%
% where A e B are borders and gamma the angle between them.
%
% Then generalized Pitagora Theorem is applied to discover if gamma is < or
% > than 90° (the sine loses the information)
%
% Output & Input of this function are:
%
% T: the considered triangle
% OppositeB: the border opposite to the angle we want extimate
% AdjacentB1: a border that is next to the angle
% AdjacentB2: the other border that is next to the angle
%
% Angle: the angle (in degrees) that we want to know
%

global TInfo

% Calculate lenghts
lAdjacentB1 = find_length_of_memorized_border (AdjacentB1);
lAdjacentB2 = find_length_of_memorized_border (AdjacentB2);
lOppositeB = find_length_of_memorized_border (OppositeB);


% Calculate sinAngle
SinAngle = 2*TInfo(T).Area / (lAdjacentB1*lAdjacentB2);

% Calculate Angle in rad
if lOppositeB^2 > (lAdjacentB1^2 + lAdjacentB2^2)
    Angle = ( pi - asin(SinAngle) );
else
    Angle = asin(SinAngle);
end

% Convert Angle to deg:
Angle = Angle*180/pi;

return