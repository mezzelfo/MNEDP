function [InitialVertexReference] = insert_vertex_in_box(InputVertex,nInputVertex)
%
% [IntialVertexReference] = insert_vertex_in_box(InputVertex,nInputVertex)
%
% This function creates a rectangular box (=Containing Box = CB) around the input vertex. The box
% is oversized with respect to the verteces. Then the function triangulates
% the box and insert the input verteces into the box, inserting new
% verteces if any of the edge border is encroached
%
% Inputs & outputs of this functions are
% 
% InputVertex(n,1): the x-coordinate of input vertex referenced ad n
% InputVertex(n,1): the x-coordinate of input vertex referenced ad n
% nInputVertex(n,1): the number of input verteces
%
% InitialVertexReference: the i-esim component of this vector says that
%   InputVertex(i) is inserted in the triangulation as the
%   InitialVertexReference(i).esim element of V array

% declaration of global 
global TV TT B V VB TInfo BInfo BCircle
global nT nV nB
global EB

global XBoxEnlarge YBoxEnlarge
global SafeMeshing

% ----------------
% Create the CB
% ----------------

% Initialize (InitialVertexReference)
InitialVertexReference = zeros(nInputVertex,1);

% Find the CB dimensions
x_max = max(InputVertex(:,1));
y_max = max(InputVertex(:,2));
x_min = min(InputVertex(:,1));
y_min = min(InputVertex(:,2));

        % (---Update EB subdomain---)
        EB.Br(1).SplittingLimits = [x_max x_min y_max y_min];

x_length = abs(x_max-x_min);
y_length = abs(y_max-y_min);

x_max = x_max + x_length * XBoxEnlarge;
y_max = y_max + y_length * YBoxEnlarge;
x_min = x_min - x_length * XBoxEnlarge;
y_min = y_min - y_length * YBoxEnlarge;



% Updating the reference geometry (ref: FIG 01)
TV = [1, 3, 2
      1, 4, 3];
  
TT = [2,-1,-1,5,3,2,3,-1,-1
      -1,-1,1,1,4,5,-1,-1,1];
  
B = [1,4,-1,2
     1,2,1,-1
     2,3,1,-1
     3,4,-1,2
     1,3,1,2 ];
 
V = [x_min,y_max
    x_max,y_max
    x_max,y_min
    x_min,y_min];

VB(1).n = 3;
VB(1).B = [1 5 2];
VB(1).V = [4 3 2];
VB(2).n = 2;
VB(2).B = [2 3];
VB(2).V = [1 3];
VB(3).n = 3;
VB(3).B = [3 5 4];
VB(3).V = [2 1 4];
VB(4).n = 2;
VB(4).B = [4 1];
VB(4).V = [3 1];

nT = 2;
nV = 4;
nB = 5;

 % Find borders info
for i = 1:5
    [BCircle(i).Circumcenter , BCircle(i).r2] = find_bordercircle_info (V(B(i,1),:),V(B(i,2),:));
end

BInfo(1,:) = [0 1 0];
BInfo(2,:) = [0 1 0];
BInfo(3,:) = [0 1 0];
BInfo(4,:) = [0 1 0];
BInfo(5,:) = [0 0 0];

% Find triangle info
TInfo(1).Area = find_area_of_memorized_tria (1);
[TInfo(1).Circumcenter(1),TInfo(1).Circumcenter(2),TInfo(1).Circumradius] = find_circuminfo_of_memorized_tria (1);
find_B_of_memorized_triangle(1);
[TInfo(1).CG]=find_CG_of_given_3_vertex ( V(TV(1,1),:), V(TV(1,2),:), V(TV(1,3),:) );

TInfo(2).Area = find_area_of_memorized_tria (2);
[TInfo(2).Circumcenter(1),TInfo(2).Circumcenter(2),TInfo(2).Circumradius] = find_circuminfo_of_memorized_tria (2);
find_B_of_memorized_triangle(2);
[TInfo(2).CG]=find_CG_of_given_3_vertex ( V(TV(2,1),:), V(TV(2,2),:), V(TV(2,3),:) );

% Modify EB (depends on EB choice)
% --------------------------------

% Update encroachable borders
for i=1:4
    encroached_borders_add (i)
end

% ------------------------------------------
% Inserct every InputVertex in the CB
% ------------------------------------------

ControlledVertex = 1; % Says which vertex I'm trying to insert in CB

% InputVertex = my_sort_ascending(InputVertex,nInputVertex,2,1);

while ControlledVertex <= nInputVertex
        
    % Research of the triangle(s) that contains the vertex
        % The search is from the bottom because InputVertex have been
        % sorted
    
    for CurrentTriangle = nT : - 1 : 1
        
        PositionVertexTriangle = check_if_in_triangle(InputVertex(ControlledVertex,1),InputVertex(ControlledVertex,2),CurrentTriangle);
        
        if PositionVertexTriangle == 0 
        
            % In this case the point belongs to the interior of triangle T
            
            EncroachedBorder = check_encroaching_border ( InputVertex(ControlledVertex,1) , InputVertex(ControlledVertex,2) ); %gives reference, -1 for none
            if EncroachedBorder == -1
                
                % In this case no border is encroached and InputVertex is accepted
                Xinsert = InputVertex(ControlledVertex,1);
                Yinsert = InputVertex(ControlledVertex,2);
                insert_vertex_given_triangle (Xinsert,Yinsert,CurrentTriangle);
                
                InitialVertexReference(ControlledVertex) = nV;
                ControlledVertex = ControlledVertex + 1;
                break
            else
                
                % In this case a border is encroached and InputVertex is rejected
                Xinsert = (V (B(EncroachedBorder,1) , 1) + V (B(EncroachedBorder,2) , 1))/2;
                Yinsert = (V (B(EncroachedBorder,1) , 2) + V (B(EncroachedBorder,2) , 2))/2;
                insert_vertex_given_border (Xinsert,Yinsert,EncroachedBorder);

                % Eventually checks for recursive encroaching
                if SafeMeshing == true
                    checks_for_recursive_encroaching (nV,1);
                end
                
                break
            end
            
        elseif PositionVertexTriangle > 0
            
            % In this case the point belongs to border PositionVertexTriangle
                        
            EncroachedBorder = check_encroaching_border ( InputVertex(ControlledVertex,1) , InputVertex(ControlledVertex,2) ); %gives reference, -1 for none
            if EncroachedBorder == -1

                % In this case no border is encroached and InputVertex is accepted
                Xinsert = InputVertex(ControlledVertex,1);
                Yinsert = InputVertex(ControlledVertex,2);
                insert_vertex_given_border (Xinsert,Yinsert,PositionVertexTriangle);
                
                InitialVertexReference(ControlledVertex) = nV;
                ControlledVertex = ControlledVertex + 1;
                break
            else
                
                % In this case a border is encroached and InputVertex is rejected
                Xinsert = (V (B(EncroachedBorder,1) , 1) + V (B(EncroachedBorder,2) , 1))/2;
                Yinsert = (V (B(EncroachedBorder,1) , 2) + V (B(EncroachedBorder,2) , 2))/2;
                insert_vertex_given_border (Xinsert,Yinsert,EncroachedBorder);
                
                % Eventually checks for recursive encroaching
                if SafeMeshing == true
                    checks_for_recursive_encroaching (nV,1);
                end

                break
            end
            
        end
        
    end % CurrentTriangle   
        
end % while --- ControlledVertex


return