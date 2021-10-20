function [] = insert_vertex_given_border (Vx,Vy,Border);
%
% [] = insert_vertex_given_triangle (Vx,Vy,Border)
%
% This function locally modify an existing grid inserting a vertex whose
% coordinates are (Vx,Vy). The vertex must belong to a border whose
% reference in the grid is Boredr (B must be known!). In a first moment the
% functions only deletes B and ajdiacent triangles and recreate two or four
% new triangles (these cases will be explained), then apply the
% myflip (*see attachments) algorithm to transform the grid into a Delaunay
% Triangolation (*see attachments)
%
% Inputs of this function are:
%
% [Vx,Vy]: the coordinates of the vertex that must be inserted;
% B: the reference to the Border
%

global TV TT VT B V VB TInfo BInfo BCircle
global nT nV nB

% Check if the vertex belongs to the boundary of the domain or not. If it
% belongs two new triangles will be creates, otherwise four

if BInfo(Border,2) == 1; %in this case the border is a boundary border
    
    % Find the reference to T
    if B(Border,3) == -1
        T = B(Border,4);
    else
        T = B(Border,3);
    end
    
    % Find references to Neighbour Triangles
    if TT(T,4) == Border
        ReferenceOfOuterRegionInTT = 1; %Says that TT(T,1) = -1
        Neighbour1 = TT(T,2);
        Neighbour2 = TT(T,3);
        Border1 = TT(T,5);
        Border2 = TT(T,6);
        ReferenceInNeighbour1 = TT(T,8);
        ReferenceInNeighbour2 = TT(T,9);
        VertexWithNeighbour1 = B(Border1,1:2); % These vertex are shared between T & N1
        VertexWithNeighbour2 = B(Border2,1:2); % These vertex are shared between T & N2
    elseif TT(T,5) == Border
        ReferenceOfOuterRegionInTT = 2; %Says that TT(T,2) = -1
        Neighbour1 = TT(T,1);
        Neighbour2 = TT(T,3);
        Border1 = TT(T,4);
        Border2 = TT(T,6);
        ReferenceInNeighbour1 = TT(T,7);
        ReferenceInNeighbour2 = TT(T,9);
        VertexWithNeighbour1 = B(Border1,1:2); % These vertex are shared between T & N1
        VertexWithNeighbour2 = B(Border2,1:2); % These vertex are shared between T & N2
    else
        ReferenceOfOuterRegionInTT = 3; %Says that TT(T,3) = -1
        Neighbour1 = TT(T,1);
        Neighbour2 = TT(T,2);
        Border1 = TT(T,4);
        Border2 = TT(T,5);
        ReferenceInNeighbour1 = TT(T,7);
        ReferenceInNeighbour2 = TT(T,8);
        VertexWithNeighbour1 = B(Border1,1:2); % These vertex are shared between T & N1
        VertexWithNeighbour2 = B(Border2,1:2); % These vertex are shared between T & N2
    end 
    
    % Find reference to BorderInfo
    BorderInfo = BInfo(Border,:);
    
    %Update vertex list
    nV = nV + 1;
    V(nV,:) = [Vx,Vy];

    % Create new references
    NewTriangle1 = T;
    NewTriangle2 = nT+1;

    NewOuterBorder1 = Border;
    NewOuterBorder2 = nB + 1;
    NewInnerBorder = nB + 2;

    % Update TV
    TV(NewTriangle1,:) = [nV,VertexWithNeighbour1];
    TV(NewTriangle2,:) = [nV,VertexWithNeighbour2];

    % Update TT in NewTriangles
    TT(NewTriangle1,:) = [Neighbour1,NewTriangle2,-1,Border1,NewInnerBorder,NewOuterBorder1,ReferenceInNeighbour1,2,-1];
    TT(NewTriangle2,:) = [Neighbour2,NewTriangle1,-1,Border2,NewInnerBorder,NewOuterBorder2,ReferenceInNeighbour2,2,-1];

    % Eventually update TT in Neighbours
    if Neighbour1 ~= -1
        TT(Neighbour1,ReferenceInNeighbour1) = NewTriangle1;
        TT(Neighbour1,ReferenceInNeighbour1+6) = 1;
    end
    if Neighbour2 ~= -1
        TT(Neighbour2,ReferenceInNeighbour2) = NewTriangle2;
        TT(Neighbour2,ReferenceInNeighbour2+6) = 1;
    end

    % Update exixsting references to B
    B(Border1,3:4)=[Neighbour1,NewTriangle1];
    B(Border2,3:4)=[Neighbour2,NewTriangle2];
    
    % Find which verteces belong to NewBorders
    if VertexWithNeighbour1(1) == VertexWithNeighbour2(1)
        CommonVertex = VertexWithNeighbour1(1);
        VertexOfOnlyNeighbour1 = VertexWithNeighbour1(2);
        VertexOfOnlyNeighbour2 = VertexWithNeighbour2(2);
    elseif VertexWithNeighbour1(1) == VertexWithNeighbour2(2)
        CommonVertex = VertexWithNeighbour1(1);
        VertexOfOnlyNeighbour1 = VertexWithNeighbour1(2);
        VertexOfOnlyNeighbour2 = VertexWithNeighbour2(1);
    elseif VertexWithNeighbour1(2) == VertexWithNeighbour2(1)
        CommonVertex = VertexWithNeighbour1(2);
        VertexOfOnlyNeighbour1 = VertexWithNeighbour1(1);
        VertexOfOnlyNeighbour2 = VertexWithNeighbour2(2);
    else
        CommonVertex = VertexWithNeighbour1(2);
        VertexOfOnlyNeighbour1 = VertexWithNeighbour1(1);
        VertexOfOnlyNeighbour2 = VertexWithNeighbour2(1);
    end

    % Creating NewBorder references
    B(NewInnerBorder,:) = [nV,CommonVertex,NewTriangle1,NewTriangle2];
    B(NewOuterBorder1,:) = [VertexOfOnlyNeighbour1,nV,NewTriangle1,-1];
    B(NewOuterBorder2,:) = [VertexOfOnlyNeighbour2,nV,NewTriangle2,-1];
    
    % Create new borders info
    [BCircle(NewInnerBorder).Circumcenter , BCircle(NewInnerBorder).r2] = find_bordercircle_info (V(B(NewInnerBorder,1),:),V(B(NewInnerBorder,2),:));
    [BCircle(NewOuterBorder1).Circumcenter , BCircle(NewOuterBorder1).r2] = find_bordercircle_info (V(B(NewOuterBorder1,1),:),V(B(NewOuterBorder1,2),:));
    [BCircle(NewOuterBorder2).Circumcenter , BCircle(NewOuterBorder2).r2] = find_bordercircle_info (V(B(NewOuterBorder2,1),:),V(B(NewOuterBorder2,2),:));
    
    % Set new borders conditions; these will eventually change
    BInfo(NewInnerBorder,:)=[0 0 0];
    BInfo(NewOuterBorder1,:) = BorderInfo;
    BInfo(NewOuterBorder2,:) = BorderInfo;
    
    % Eventually modify encroachable list
    if BorderInfo(2) == 1
        encroached_borders_add (NewOuterBorder2)
    end
        
    %Update references to VB
    VB(CommonVertex).n = VB(CommonVertex).n + 1;
    VB(CommonVertex).V(VB(CommonVertex).n) = nV;
    VB(CommonVertex).B(VB(CommonVertex).n) = NewInnerBorder;
    
    for i = 1: VB(VertexOfOnlyNeighbour1).n 
        if VB(VertexOfOnlyNeighbour1).V(i) == VertexOfOnlyNeighbour2
            VB(VertexOfOnlyNeighbour1).V(i) = nV;
            VB(VertexOfOnlyNeighbour1).B(i) = NewOuterBorder1;
            break
        end
    end
    
    for i = 1: VB(VertexOfOnlyNeighbour2).n 
        if VB(VertexOfOnlyNeighbour2).V(i) == VertexOfOnlyNeighbour1
            VB(VertexOfOnlyNeighbour2).V(i) = nV;
            VB(VertexOfOnlyNeighbour2).B(i) = NewOuterBorder2;
            break
        end
    end

    VB(nV).n = 3;
    VB(nV).V = [VertexOfOnlyNeighbour1,VertexOfOnlyNeighbour2,CommonVertex];
    VB(nV).B = [NewOuterBorder1,NewOuterBorder2,NewInnerBorder];

    % Find info about NewTriangles
    TInfo(NewTriangle1).Area = find_area_of_memorized_tria (NewTriangle1);
    [TInfo(NewTriangle1).Circumcenter(1),TInfo(NewTriangle1).Circumcenter(2),TInfo(NewTriangle1).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle1);
    find_B_of_memorized_triangle(NewTriangle1);
    [TInfo(NewTriangle1).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle1,1),:),V(TV(NewTriangle1,2),:),V(TV(NewTriangle1,3),:) );

    
    TInfo(NewTriangle2).Area = find_area_of_memorized_tria (NewTriangle2);
    [TInfo(NewTriangle2).Circumcenter(1),TInfo(NewTriangle2).Circumcenter(2),TInfo(NewTriangle2).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle2);
    find_B_of_memorized_triangle(NewTriangle1);
    [TInfo(NewTriangle2).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle2,1),:),V(TV(NewTriangle2,2),:),V(TV(NewTriangle2,3),:) );

    
    %Update totals
    nT = nT + 1;
    nB = nB + 2;
    
    % Check if some of the new borders are kidnapped
    check_if_kidnapped_are_found(NewInnerBorder);
    check_if_kidnapped_are_found(NewOuterBorder1);
    check_if_kidnapped_are_found(NewOuterBorder2);
    
    % Myflip actual grid
    myflip (NewTriangle1,Neighbour1,Border1);
    myflip (NewTriangle2,Neighbour2,Border2);
    
else %in this case the border is an inner border
    
    % Find the reference to T
    T1 =B(Border,3);
    T2 =B(Border,4);
    
    % Find references to Neighbour Triangles 11 e 12
    if TT(T1,1) == T2
        Neighbour11 = TT(T1,2);
        Neighbour12 = TT(T1,3);
        Border11 = TT(T1,5);
        Border12 = TT(T1,6);
        ReferenceInNeighbour11 = TT(T1,8);
        ReferenceInNeighbour12 = TT(T1,9);
        VertexWithNeighbour11 = B(Border11,1:2); % These vertex are shared between T & N11
        VertexWithNeighbour12 = B(Border12,1:2); % These vertex are shared between T & N12
        ReferenceOfT1InT2 = TT(T1,7); %Says that TT(T2,ReferenceOfT1InT2) = T1
    elseif TT(T1,2) == T2
        Neighbour11 = TT(T1,1);
        Neighbour12 = TT(T1,3);
        Border11 = TT(T1,4);
        Border12 = TT(T1,6);
        ReferenceInNeighbour11 = TT(T1,7);
        ReferenceInNeighbour12 = TT(T1,9);
        VertexWithNeighbour11 = B(Border11,1:2); % These vertex are shared between T & N11
        VertexWithNeighbour12 = B(Border12,1:2); % These vertex are shared between T & N12
        ReferenceOfT1InT2 = TT(T1,8); %Says that TT(T2,ReferenceOfT1InT2) = T1
    else
        Neighbour11 = TT(T1,1);
        Neighbour12 = TT(T1,2);
        Border11 = TT(T1,4);
        Border12 = TT(T1,5);
        ReferenceInNeighbour11 = TT(T1,7);
        ReferenceInNeighbour12 = TT(T1,8);
        VertexWithNeighbour11 = B(Border11,1:2); % These vertex are shared between T & N11
        VertexWithNeighbour12 = B(Border12,1:2); % These vertex are shared between T & N12
        ReferenceOfT1InT2 = TT(T1,9); %Says that TT(T2,ReferenceOfT1InT2) = T1
    end 

    % Find which verteces belong to new borders
    if VertexWithNeighbour11(1) == VertexWithNeighbour12(1)
        FarVertex1 = VertexWithNeighbour11(1);
        NearVertex14 = VertexWithNeighbour11(2);
        NearVertex23 = VertexWithNeighbour12(2);
    elseif VertexWithNeighbour11(1) == VertexWithNeighbour12(2)
        FarVertex1 = VertexWithNeighbour11(1);
        NearVertex14 = VertexWithNeighbour11(2);
        NearVertex23 = VertexWithNeighbour12(1);
    elseif VertexWithNeighbour11(2) == VertexWithNeighbour12(1)
        FarVertex1 = VertexWithNeighbour11(2);
        NearVertex14 = VertexWithNeighbour11(1);
        NearVertex23 = VertexWithNeighbour12(2);
    else
        FarVertex1 = VertexWithNeighbour11(2);
        NearVertex14 = VertexWithNeighbour11(1);
        NearVertex23 = VertexWithNeighbour12(1);
    end

    
    % Find references to Neighbour Triangles 23 e 24
    
    if ReferenceOfT1InT2 == 1 %Find reference in T2
        OtherReference1 = 2;
        OtherReference2 = 3;
    elseif ReferenceOfT1InT2 == 2
        OtherReference1 = 1;
        OtherReference2 = 3;
    else
        OtherReference1 = 1;
        OtherReference2 = 2;
    end
    
    if B( TT(T2,OtherReference1 + 3) , 1 ) == NearVertex23 | B( TT(T2,OtherReference1 + 3) , 2 ) == NearVertex23
        Neighbour23 = TT(T2,OtherReference1);
        Neighbour24 = TT(T2,OtherReference2);
        Border23 = TT(T2,OtherReference1 + 3);
        Border24 = TT(T2,OtherReference2 + 3);
        ReferenceInNeighbour23 = TT(T2,OtherReference1 + 6);
        ReferenceInNeighbour24 = TT(T2,OtherReference2 + 6);
        VertexWithNeighbour23 = B(Border23,1:2); % These vertex are shared between T & N11
        VertexWithNeighbour24 = B(Border24,1:2); % These vertex are shared between T & N12
    else
        Neighbour23 = TT(T2,OtherReference2);
        Neighbour24 = TT(T2,OtherReference1);
        Border23 = TT(T2,OtherReference2 + 3);
        Border24 = TT(T2,OtherReference1 + 3);
        ReferenceInNeighbour23 = TT(T2,OtherReference2 + 6);
        ReferenceInNeighbour24 = TT(T2,OtherReference1 + 6);
        VertexWithNeighbour23 = B(Border23,1:2); % These vertex are shared between T & N11
        VertexWithNeighbour24 = B(Border24,1:2); % These vertex are shared between T & N12
    end
    
    % Find FarVertex2
    if NearVertex14 == VertexWithNeighbour24(1)
        FarVertex2 = VertexWithNeighbour24(2);
    else
        FarVertex2 = VertexWithNeighbour24(1);
    end
    
    % Find reference to BorderInfo
    BorderInfo = BInfo(Border,:);
    
    %Update vertex list
    nV = nV + 1;
    V(nV,:) = [Vx,Vy];

    % Create new references
    NewTriangle11 = T1;
    NewTriangle12 = T2;
    NewTriangle23 = nT + 1;
    NewTriangle24 = nT + 2;
    
    NewFarBorder1 = Border;
    NewFarBorder2 = nB + 1;
    NewNearBorder14 = nB + 2;
    NewNearBorder23 = nB + 3;

    % Update TV
    TV(NewTriangle11,:) = [nV,VertexWithNeighbour11];
    TV(NewTriangle12,:) = [nV,VertexWithNeighbour12];
    TV(NewTriangle23,:) = [nV,VertexWithNeighbour23];
    TV(NewTriangle24,:) = [nV,VertexWithNeighbour24];

    % Update TT in NewTriangles
    TT(NewTriangle11,:) = [Neighbour11,NewTriangle24,NewTriangle12,Border11,NewNearBorder14,NewFarBorder1,ReferenceInNeighbour11,3,2];
    TT(NewTriangle12,:) = [Neighbour12,NewTriangle11,NewTriangle23,Border12,NewFarBorder1,NewNearBorder23,ReferenceInNeighbour12,3,2];
    TT(NewTriangle23,:) = [Neighbour23,NewTriangle12,NewTriangle24,Border23,NewNearBorder23,NewFarBorder2,ReferenceInNeighbour23,3,2];
    TT(NewTriangle24,:) = [Neighbour24,NewTriangle23,NewTriangle11,Border24,NewFarBorder2,NewNearBorder14,ReferenceInNeighbour24,3,2];
    
    
    % Eventually update TT in Neighbours
    if Neighbour11 ~= -1
        TT(Neighbour11,ReferenceInNeighbour11) = NewTriangle11;
        TT(Neighbour11,ReferenceInNeighbour11+6) = 1;
    end
    if Neighbour12 ~= -1
        TT(Neighbour12,ReferenceInNeighbour12) = NewTriangle12;
        TT(Neighbour12,ReferenceInNeighbour12+6) = 1;
    end
    if Neighbour23 ~= -1
        TT(Neighbour23,ReferenceInNeighbour23) = NewTriangle23;
        TT(Neighbour23,ReferenceInNeighbour23+6) = 1;
    end
    if Neighbour24 ~= -1
        TT(Neighbour24,ReferenceInNeighbour24) = NewTriangle24;
        TT(Neighbour24,ReferenceInNeighbour24+6) = 1;
    end

    % Update existing references to B
    B(Border11,3:4)=[Neighbour11,NewTriangle11];
    B(Border12,3:4)=[Neighbour12,NewTriangle12];
    B(Border23,3:4)=[Neighbour23,NewTriangle23];
    B(Border24,3:4)=[Neighbour24,NewTriangle24];

   
    % Creating new borders references
    B(NewNearBorder14,:) = [nV,NearVertex14,NewTriangle11,NewTriangle24];
    B(NewFarBorder1,:) = [nV,FarVertex1,NewTriangle11,NewTriangle12];
    B(NewNearBorder23,:) = [nV,NearVertex23,NewTriangle12,NewTriangle23];
    B(NewFarBorder2,:) = [nV,FarVertex2,NewTriangle23,NewTriangle24];

    % Create new borders info
    [BCircle(NewNearBorder14).Circumcenter , BCircle(NewNearBorder14).r2] = find_bordercircle_info (V(B(NewNearBorder14,1),:),V(B(NewNearBorder14,2),:));
    [BCircle(NewNearBorder23).Circumcenter , BCircle(NewNearBorder23).r2] = find_bordercircle_info (V(B(NewNearBorder23,1),:),V(B(NewNearBorder23,2),:));
    [BCircle(NewFarBorder1).Circumcenter , BCircle(NewFarBorder1).r2] = find_bordercircle_info (V(B(NewFarBorder1,1),:),V(B(NewFarBorder1,2),:));
    [BCircle(NewFarBorder2).Circumcenter , BCircle(NewFarBorder2).r2] = find_bordercircle_info (V(B(NewFarBorder2,1),:),V(B(NewFarBorder2,2),:));
    
    % Set new borders conditions; these will eventually change
    BInfo(NewFarBorder1,:) = [0 0 0];
    BInfo(NewFarBorder2,:) = [0 0 0];
    BInfo(NewNearBorder14,:) = BorderInfo;
    BInfo(NewNearBorder23,:) = BorderInfo;
    
    %Update references to VB
    VB(FarVertex1).n = VB(FarVertex1).n + 1;
    VB(FarVertex1).V(VB(FarVertex1).n) = nV;
    VB(FarVertex1).B(VB(FarVertex1).n) = NewFarBorder1;
    
    VB(FarVertex2).n = VB(FarVertex2).n + 1;
    VB(FarVertex2).V(VB(FarVertex2).n) = nV;
    VB(FarVertex2).B(VB(FarVertex2).n) = NewFarBorder2;

    for i = 1: VB(NearVertex23).n 
        if VB(NearVertex23).V(i) == NearVertex14
            VB(NearVertex23).V(i) = nV;
            VB(NearVertex23).B(i) = NewNearBorder23;
            break
        end
    end
    
    for i = 1: VB(NearVertex14).n 
        if VB(NearVertex14).V(i) == NearVertex23
            VB(NearVertex14).V(i) = nV;
            VB(NearVertex14).B(i) = NewNearBorder14;
            break
        end
    end

    VB(nV).n = 4;
    VB(nV).V = [NearVertex14,FarVertex1,NearVertex23,FarVertex2];
    VB(nV).B = [NewNearBorder14,NewFarBorder1,NewNearBorder23,NewFarBorder2];

    % Find info about NewTriangles
    TInfo(NewTriangle11).Area = find_area_of_memorized_tria (NewTriangle11);
    [TInfo(NewTriangle11).Circumcenter(1),TInfo(NewTriangle11).Circumcenter(2),TInfo(NewTriangle11).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle11);
    find_B_of_memorized_triangle(NewTriangle11);
    [TInfo(NewTriangle11).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle11,1),:),V(TV(NewTriangle11,2),:),V(TV(NewTriangle11,3),:) );

    
    TInfo(NewTriangle12).Area = find_area_of_memorized_tria (NewTriangle12);
    [TInfo(NewTriangle12).Circumcenter(1),TInfo(NewTriangle12).Circumcenter(2),TInfo(NewTriangle12).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle12);
    find_B_of_memorized_triangle(NewTriangle12);
    [TInfo(NewTriangle12).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle12,1),:),V(TV(NewTriangle12,2),:),V(TV(NewTriangle12,3),:) );
    
    TInfo(NewTriangle23).Area = find_area_of_memorized_tria (NewTriangle23);
    [TInfo(NewTriangle23).Circumcenter(1),TInfo(NewTriangle23).Circumcenter(2),TInfo(NewTriangle23).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle23);
    find_B_of_memorized_triangle(NewTriangle23);
    [TInfo(NewTriangle23).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle23,1),:),V(TV(NewTriangle23,2),:),V(TV(NewTriangle23,3),:) );
    
    TInfo(NewTriangle24).Area = find_area_of_memorized_tria (NewTriangle24);
    [TInfo(NewTriangle24).Circumcenter(1),TInfo(NewTriangle24).Circumcenter(2),TInfo(NewTriangle24).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle24);
    find_B_of_memorized_triangle(NewTriangle24);
    [TInfo(NewTriangle24).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle24,1),:),V(TV(NewTriangle24,2),:),V(TV(NewTriangle24,3),:) );
    
    %Update totals
    nT = nT + 2;
    nB = nB + 3;
    
    % Check if some of the new borders are kidnapped
    check_if_kidnapped_are_found(NewNearBorder14);
    check_if_kidnapped_are_found(NewFarBorder1);
    check_if_kidnapped_are_found(NewNearBorder23);
    check_if_kidnapped_are_found(NewFarBorder2);

    % Myflip actual grid
    myflip (NewTriangle11,Neighbour11,Border11);
    myflip (NewTriangle12,Neighbour12,Border12);
    myflip (NewTriangle23,Neighbour23,Border23);
    myflip (NewTriangle24,Neighbour24,Border24);

    
    
end


return
