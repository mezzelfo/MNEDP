function [flag] = check_if_in_triangle (xP,yP,T);
%
% [flag] = check_if_in_triangle (x,y,T)
%
% This function controls if the given point belongs to referenced triangle,
% using the shape functions:
%
% A=[x1,x2,x3
%    y1,y2,y3
%    1,1,1]
% v=[xP
%    yP
%    1]
% lambda=A\v;
%
% Output & inputs of this function are:
%
% flag: a boolean that may assume the value:
%       0 if the point belongs to the interior of the triangle
%       n if the point belongs to the boundary of the triangle, where n is
%         the reference to the belonging border 
%       -1 if the point lyes outside the triangle
%
% xP: the x-coordinate of the point to be checked;
% yC: the y-coordinate of the point to be checked;
% T: the reference to the chosen triangle
%

global TV V

x1 = V(TV(T,1),1);
x2 = V(TV(T,2),1);
x3 = V(TV(T,3),1);
y1 = V(TV(T,1),2);
y2 = V(TV(T,2),2);
y3 = V(TV(T,3),2);

% % A=[x1,x2,x3
% %    y1,y2,y3
% %    1,1,1];
v=[xP
   yP
   1];
lambda=[y2-y3,-x2+x3,x2*y3-x3*y2 ; -y1+y3,x1-x3,-x1*y3+x3*y1 ; y1-y2,-x1+x2,x1*y2-x2*y1]*v/(x1*y2-x1*y3-x2*y1+x2*y3+x3*y1-x3*y2);

flag = -1;
% % % if ((0<=lambda(1) & lambda(1)<=1) & (0<=lambda(2) & lambda(2)<=1) & (0<=lambda(3) & lambda(3)<=1))
% % %     if (my_equal(lambda(1),0)==1)
% % %         flag = TT(T,4);
% % %     elseif(my_equal(lambda(2),0)==1)
% % %         flag = TT(T,5);
% % %     elseif (my_equal(lambda(3),0)==1)
% % %         flag = TT(T,6);
% % %     else
% % %         flag = 0;
% % %     end
% % % end

if ((0<=lambda(1) & lambda(1)<=1) & (0<=lambda(2) & lambda(2)<=1) & (0<=lambda(3) & lambda(3)<=1))
    flag = 0;
    if  my_equal(lambda(1),0) | my_equal(lambda(2),0) | my_equal(lambda(3),0)
        
        global B TT
        for I= 4:6
            Position = find_alignment_of_given_vertex( B(TT(T,I),1) , B(TT(T,I),2) , [xP yP] );
            if Position == 0
                flag = TT(T,I);
                break
            end
        end
        
    end
    
end


return