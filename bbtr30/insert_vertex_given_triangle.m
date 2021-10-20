function [] = insert_vertex_given_triangle (Vx,Vy,T);
%
% [] = insert_vertex_given_triangle (Vx,Vy,T)
%
% This function locally modify an existing grid inserting a vertex whose
% coordinates are (Vx,Vy). The vertex must lye inside a triangle whose
% reference in the grid is T (T must be known!). In a first moment the
% functions only deletes T and recreate three new triangles, then apply the
% myflip (*see attachments) algorithm to transform the grid into a Delaunay
% Triangolation (*see attachments)
%
% Inputs of this function are:
%
% [Vx,Vy]: the coordinates of the vertex that must be inserted;
% T: the reference to the triangle
%

global TV TT VT B V VB TInfo BInfo BCircle
global nT nV nB

% Make old references explicit
Neighbour1 = TT(T,1);
Neighbour2 = TT(T,2);
Neighbour3 = TT(T,3);
Border1 = TT(T,4);
Border2 = TT(T,5);
Border3 = TT(T,6);
ReferenceInNeighbour1 = TT(T,7);
ReferenceInNeighbour2 = TT(T,8);
ReferenceInNeighbour3 = TT(T,9);
VertexWithNeighbour1 = B(Border1,1:2); % These vertex are shared between T & N1
VertexWithNeighbour2 = B(Border2,1:2); % These vertex are shared between T & N2
VertexWithNeighbour3 = B(Border3,1:2); % These vertex are shared between T & N3


%Update vertex list
nV = nV + 1;
V(nV,:) = [Vx,Vy];

% Create new references
NewTriangle1 = T;
NewTriangle2 = nT+1;
NewTriangle3 = nT+2;

NewBorder1 = nB + 1;
NewBorder2 = nB + 2;
NewBorder3 = nB + 3;

% Update TV
TV(NewTriangle1,:) = [nV,VertexWithNeighbour1];
TV(NewTriangle2,:) = [nV,VertexWithNeighbour2];
TV(NewTriangle3,:) = [nV,VertexWithNeighbour3];

% Update TT in NewTriangles
TT(NewTriangle1,:) = [NewTriangle2,NewTriangle3,Neighbour1,NewBorder1,NewBorder3,Border1,2,1,ReferenceInNeighbour1];
TT(NewTriangle2,:) = [NewTriangle3,NewTriangle1,Neighbour2,NewBorder2,NewBorder1,Border2,2,1,ReferenceInNeighbour2];
TT(NewTriangle3,:) = [NewTriangle1,NewTriangle2,Neighbour3,NewBorder3,NewBorder2,Border3,2,1,ReferenceInNeighbour3];

% Eventually update TT in Neighbours
if Neighbour1 ~= -1
    TT(Neighbour1,ReferenceInNeighbour1) = NewTriangle1;
    TT(Neighbour1,ReferenceInNeighbour1+6) = 3;
end
if Neighbour2 ~= -1
    TT(Neighbour2,ReferenceInNeighbour2) = NewTriangle2;
    TT(Neighbour2,ReferenceInNeighbour2+6) = 3;
end
if Neighbour3 ~= -1
    TT(Neighbour3,ReferenceInNeighbour3) = NewTriangle3;
    TT(Neighbour3,ReferenceInNeighbour3+6) = 3;
end

% Update exixsting references to B
B(Border1,3:4)=[Neighbour1,NewTriangle1];
B(Border2,3:4)=[Neighbour2,NewTriangle2];
B(Border3,3:4)=[Neighbour3,NewTriangle3];

% Find which verteces belong to NewBorders
if VertexWithNeighbour1(1) == VertexWithNeighbour2(1)
    VertexOfNewBorder1 = VertexWithNeighbour1(1);
    VertexOfNewBorder2 = VertexWithNeighbour2(2);
    VertexOfNewBorder3 = VertexWithNeighbour1(2);
elseif VertexWithNeighbour1(1) == VertexWithNeighbour2(2)
    VertexOfNewBorder1 = VertexWithNeighbour1(1);
    VertexOfNewBorder2 = VertexWithNeighbour2(1);
    VertexOfNewBorder3 = VertexWithNeighbour1(2);
elseif VertexWithNeighbour1(2) == VertexWithNeighbour2(1)
    VertexOfNewBorder1 = VertexWithNeighbour1(2);
    VertexOfNewBorder2 = VertexWithNeighbour2(2);
    VertexOfNewBorder3 = VertexWithNeighbour1(1);
else
    VertexOfNewBorder1 = VertexWithNeighbour1(2);
    VertexOfNewBorder2 = VertexWithNeighbour2(1);
    VertexOfNewBorder3 = VertexWithNeighbour1(1);
end

% Creating NewBorder references
B(NewBorder1,:) = [nV,VertexOfNewBorder1,NewTriangle1,NewTriangle2];
B(NewBorder2,:) = [nV,VertexOfNewBorder2,NewTriangle2,NewTriangle3];
B(NewBorder3,:) = [nV,VertexOfNewBorder3,NewTriangle1,NewTriangle3];

% Set NewBorder conditions
BInfo(NewBorder1,:)=[0 0 0];
BInfo(NewBorder2,:)=[0 0 0];
BInfo(NewBorder3,:)=[0 0 0];

% Create NewBorder info
[BCircle(NewBorder1).Circumcenter , BCircle(NewBorder1).r2] = find_bordercircle_info (V(B(NewBorder1,1),:),V(B(NewBorder1,2),:));
[BCircle(NewBorder2).Circumcenter , BCircle(NewBorder2).r2] = find_bordercircle_info (V(B(NewBorder2,1),:),V(B(NewBorder2,2),:));
[BCircle(NewBorder3).Circumcenter , BCircle(NewBorder3).r2] = find_bordercircle_info (V(B(NewBorder3,1),:),V(B(NewBorder3,2),:));

%Update references to VB
VB(VertexOfNewBorder1).n = VB(VertexOfNewBorder1).n + 1;
VB(VertexOfNewBorder1).V(VB(VertexOfNewBorder1).n) = nV;
VB(VertexOfNewBorder1).B(VB(VertexOfNewBorder1).n) = NewBorder1;

VB(VertexOfNewBorder2).n = VB(VertexOfNewBorder2).n + 1;
VB(VertexOfNewBorder2).V(VB(VertexOfNewBorder2).n) = nV;
VB(VertexOfNewBorder2).B(VB(VertexOfNewBorder2).n) = NewBorder2;

VB(VertexOfNewBorder3).n = VB(VertexOfNewBorder3).n + 1;
VB(VertexOfNewBorder3).V(VB(VertexOfNewBorder3).n) = nV;
VB(VertexOfNewBorder3).B(VB(VertexOfNewBorder3).n) = NewBorder3;

VB(nV).n = 3;
VB(nV).V = [VertexOfNewBorder1,VertexOfNewBorder2,VertexOfNewBorder3];
VB(nV).B = [NewBorder1,NewBorder2,NewBorder3];

% Find info about NewTriangles
TInfo(NewTriangle1).Area = find_area_of_memorized_tria (NewTriangle1);
[TInfo(NewTriangle1).Circumcenter(1),TInfo(NewTriangle1).Circumcenter(2),TInfo(NewTriangle1).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle1);
find_B_of_memorized_triangle(NewTriangle1);
[TInfo(NewTriangle1).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle1,1),:),V(TV(NewTriangle1,2),:),V(TV(NewTriangle1,3),:) );


TInfo(NewTriangle2).Area = find_area_of_memorized_tria (NewTriangle2);
[TInfo(NewTriangle2).Circumcenter(1),TInfo(NewTriangle2).Circumcenter(2),TInfo(NewTriangle2).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle2);
find_B_of_memorized_triangle(NewTriangle2);
[TInfo(NewTriangle2).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle2,1),:),V(TV(NewTriangle2,2),:),V(TV(NewTriangle2,3),:) );


TInfo(NewTriangle3).Area = find_area_of_memorized_tria (NewTriangle3);
[TInfo(NewTriangle3).Circumcenter(1),TInfo(NewTriangle3).Circumcenter(2),TInfo(NewTriangle3).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle3);
find_B_of_memorized_triangle(NewTriangle3);
[TInfo(NewTriangle3).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle3,1),:),V(TV(NewTriangle3,2),:),V(TV(NewTriangle3,3),:) );


%Update totals
nT = nT + 2;
nB = nB + 3;

% Check if some of the new borders are kidnapped
check_if_kidnapped_are_found(NewBorder1);
check_if_kidnapped_are_found(NewBorder2);
check_if_kidnapped_are_found(NewBorder3);

% Myflip actual grid
myflip (NewTriangle1,Neighbour1,Border1);
myflip (NewTriangle2,Neighbour2,Border2);
myflip (NewTriangle3,Neighbour3,Border3);

return
