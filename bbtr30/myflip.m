function [] = myflip (T1,T2,Border);
%
% [] = myflip (T1,T2,Border)
%
% This function locally modify an existing grid, checking if a border must
% be myflipped (*see attachments) to recreate a Delaunay triangulation (*see
% attachments). The function controls if a particular vertex, FarVertex2
% (see figure 5), lyes within border circumcircle (*see attachments), and
% in positive case myflips the border. Then the function checks if newborn
% triangles must be myflipped again (the functionis recursive)
%
% Inputs of this function are:
%
% T1: the triangle that must be controlled
% T2: this is the other triangle that may be myflipped. "The myflip wawe" goes
%     in this direction
% Border: the reference to the Border that may be myflipped

if T2 ~= -1
    
    global B TV V TInfo
    
    % Determine FarVertex2
    if B(Border,1) == TV(T2,1)
        if B(Border,2) == TV(T2,2)
            FarVertex2 = TV(T2,3);
        else
            FarVertex2 = TV(T2,2);
        end
    elseif B(Border,2) == TV(T2,1)
        if B(Border,1) == TV(T2,2)
            FarVertex2 = TV(T2,3);
        else
            FarVertex2 = TV(T2,2);
        end
    else
        FarVertex2 = TV(T2,1);
    end
    
    % Discover if VarVertex2 lyes within the bordercircle of border
    MyflipPosition = check_if_in_circle (V(FarVertex2,1),V(FarVertex2,2),TInfo(T1).Circumcenter,TInfo(T1).Circumradius);
    
    if MyflipPosition <= 0 %Myflip borders
        
        % Declare other global variables
        
        global TT VB BCircle BInfo
        
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
        
         % Find FarVertex1,NearVertex14 % NearVertex23
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
        
        % Find references to Neighbour Triangles 23 e 24
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
        
       
        % Create new references
        NewTriangle1 = T1;
        NewTriangle2 = T2;
        % Border is not renamed
        
        % Update TV
        TV(NewTriangle1,:) = [FarVertex1,NearVertex14,FarVertex2];
        TV(NewTriangle2,:) = [FarVertex1,FarVertex2,NearVertex23];
        
        % Update TT in NewTriangles
        TT(NewTriangle1,:) = [Neighbour11,Neighbour24,NewTriangle2,Border11,Border24,Border,ReferenceInNeighbour11,ReferenceInNeighbour24,3];
        TT(NewTriangle2,:) = [Neighbour12,Neighbour23,NewTriangle1,Border12,Border23,Border,ReferenceInNeighbour12,ReferenceInNeighbour23,3];
        
        % Eventually update TT in Neighbours
        if Neighbour11 ~= -1
            TT(Neighbour11,ReferenceInNeighbour11) = NewTriangle1;
            TT(Neighbour11,ReferenceInNeighbour11+6) = 1;
        end
        if Neighbour12 ~= -1
            TT(Neighbour12,ReferenceInNeighbour12) = NewTriangle2;
            TT(Neighbour12,ReferenceInNeighbour12+6) = 1;
        end
        if Neighbour23 ~= -1
            TT(Neighbour23,ReferenceInNeighbour23) = NewTriangle2;
            TT(Neighbour23,ReferenceInNeighbour23+6) = 2;
        end
        if Neighbour24 ~= -1
            TT(Neighbour24,ReferenceInNeighbour24) = NewTriangle1;
            TT(Neighbour24,ReferenceInNeighbour24+6) = 2;
        end
        
        % Update existing references to B
        B(Border11,3:4)=[Neighbour11,NewTriangle1];
        B(Border12,3:4)=[Neighbour12,NewTriangle2];
        B(Border23,3:4)=[Neighbour23,NewTriangle2];
        B(Border24,3:4)=[Neighbour24,NewTriangle1];
        
        % Creating new borders references
        if BInfo(Border,1) ~= 0
            border_kidnap (Border,NewTriangle1);
        end
        B(Border,:) = [FarVertex1,FarVertex2,NewTriangle1,NewTriangle2];
        BInfo(Border,:) = [0 0 0];
        
        % Find NewBorder info
        [BCircle(Border).Circumcenter , BCircle(Border).r2] = find_bordercircle_info (V(B(Border,1),:),V(B(Border,2),:));
        
        %Update references to VB
        VB(FarVertex1).n = VB(FarVertex1).n + 1;
        VB(FarVertex1).V(VB(FarVertex1).n) = FarVertex2;
        VB(FarVertex1).B(VB(FarVertex1).n) = Border;
        
        VB(FarVertex2).n = VB(FarVertex2).n + 1;
        VB(FarVertex2).V(VB(FarVertex2).n) = FarVertex1;
        VB(FarVertex2).B(VB(FarVertex2).n) = Border;
        
        for i = 1: VB(NearVertex23).n 
            if VB(NearVertex23).V(i) == NearVertex14
                [VB(NearVertex23).V] = array_remove_component (VB(NearVertex23).V,VB(NearVertex23).n,i);
                [VB(NearVertex23).B] = array_remove_component (VB(NearVertex23).B,VB(NearVertex23).n,i);
                break
            end
        end
        VB(NearVertex23).n = VB(NearVertex23).n - 1;
        
        for i = 1: VB(NearVertex14).n 
            if VB(NearVertex14).V(i) == NearVertex23
                [VB(NearVertex14).V] = array_remove_component (VB(NearVertex14).V,VB(NearVertex14).n,i);
                [VB(NearVertex14).B] = array_remove_component (VB(NearVertex14).B,VB(NearVertex14).n,i);
                break
            end
        end
        VB(NearVertex14).n = VB(NearVertex14).n - 1;
        
        % Find info about NewTriangles
        TInfo(NewTriangle1).Area = find_area_of_memorized_tria (NewTriangle1);
        [TInfo(NewTriangle1).Circumcenter(1),TInfo(NewTriangle1).Circumcenter(2),TInfo(NewTriangle1).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle1);
        find_B_of_memorized_triangle(NewTriangle1);
        [TInfo(NewTriangle1).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle1,1),:),V(TV(NewTriangle1,2),:),V(TV(NewTriangle1,3),:) );
        
        TInfo(NewTriangle2).Area = find_area_of_memorized_tria (NewTriangle2);
        [TInfo(NewTriangle2).Circumcenter(1),TInfo(NewTriangle2).Circumcenter(2),TInfo(NewTriangle2).Circumradius] = find_circuminfo_of_memorized_tria (NewTriangle2);
        find_B_of_memorized_triangle(NewTriangle2)
        [TInfo(NewTriangle2).CG]=find_CG_of_given_3_vertex ( V(TV(NewTriangle2,1),:),V(TV(NewTriangle2,2),:),V(TV(NewTriangle2,3),:) );
        
        % Check if some of the new borders are kidnapped
        check_if_kidnapped_are_found(Border);
        
        % Myflip actual grid
        myflip (NewTriangle1,Neighbour24,Border24);
        myflip (NewTriangle2,Neighbour23,Border23);
        
    end

end
    
return
