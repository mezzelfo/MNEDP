function [Flag] = check_if_in_circle (Vx,Vy,Center,r2);
%
% [Flag] = check_if_in_circle (Vx,Vy,Center,r2)
%
% This function controls if the given point (Vx,Vy) )belongs to a specific circle,
% putting values into circle equation. The check is made using setted
% tolerance
%
% Inputs & outputs of this functions are:
%
% flag: a boolean that may assume the value:
%       -1 if the point lyes inside the circle
%       0 if the point belongs to the boundary of the circle 
%       1 if the point lyes outside the circle
%
% [Vx,Vy]: the coordinates of the points to be checked
% [Center(1),Center(2)]: the coordinates of the center of the circle
% r2: the squadred radius of the circle
%

global tolerance

position = (Vx-Center(1))^2 + (Vy-Center(2))^2 - r2;
if position > tolerance
    Flag = 1;
elseif position < -tolerance
    Flag = -1;
else
    Flag = 0;
end

return