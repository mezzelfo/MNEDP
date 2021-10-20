function [xC,yC,r2]=find_circuminfo_of_memorized_tria (T);
%
% [xC,yC,r2] = find_circuminfo_of_memorized_tria (T)
%
% This function considers the circumcircle of the given triangle (passed as
% reference). Used relations are:
%
% A=[x1,y1,1
%    x2,y2,1
%    x3,y3,1]
% v=[-x1^2-y1^2,-x2^2-y2^2,-x3^2-y3^2]'
% x=A\v;
% xC=-x(1)/2
% yC=-x(2)/2
% r2=xC^2+yC^2-x(3)
%
% these are developped with symbolic method
%
% Output of this function are:
%
% xC: the x-coordinate of the circumcenter;
% yC: the y-coordinate of the circumcenter;
% r2: the square of the radius of the circumcenter (circumradius)
%
% T: the reference to the chosen triangle
%

global TV V

x1 = V(TV(T,1),1);
x2 = V(TV(T,2),1);
x3 = V(TV(T,3),1);
y1 = V(TV(T,1),2);
y2 = V(TV(T,2),2);
y3 = V(TV(T,3),2);

det = x1*y2-x1*y3-x2*y1+x2*y3+x3*y1-x3*y2;
xC = -((y2-y3)*(-x1^2-y1^2)+(-y1+y3)*(-x2^2-y2^2)+(y1-y2)*(-x3^2-y3^2))/2/det;
yC = -((-x2+x3)*(-x1^2-y1^2)+(x1-x3)*(-x2^2-y2^2)+(-x1+x2)*(-x3^2-y3^2))/2/det;
r2=xC^2+yC^2-((x2*y3-x3*y2)*(-x1^2-y1^2)+(-x1*y3+x3*y1)*(-x2^2-y2^2)+(x1*y2-x2*y1)*(-x3^2-y3^2))/det;

return