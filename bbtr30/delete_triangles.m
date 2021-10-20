function [] = delete_triangles (TriangleToBeDeleted,nTriangleToBeDeleted);
%
% [] = delete_triangles (TriangleToBeDeleted,nTriangleToBeDeleted);
%
% This function deletes triangles insert in the parameter array. If it is
% necessary the function also eliminates borders and verteces and updates
% BorderInfo array.
%
% Input of this function are:
% TriangleToBeDeleted: the list of triangles that must be deleted
% nTriangleToBeDeleted: the number of triangles that must be deletetd
%

global TV TT B V VB TInfo BInfo BCircle
global nT nV nB

global EncroachableBorders nEncroachableBorders

% Define useful variables
i = nTriangleToBeDeleted;


while i > 0
    
    iT = TriangleToBeDeleted (i);
    
    % Find which triangle will fill the hole in the list
    if iT == nT
        %iT is last triangle in the list 
        TMoved = 0;
    else
        TMoved = nT;
    end
    
    % ------------------------------------------------
    % Delete reference in NeigbourTriangle and Border
    % ------------------------------------------------
    
    
    % Initialize references
    % ----------------------
    
    nBorderToDelete = 0;
    nVertexToDelete = 0;
    BorderToDelete = [];
    VertexToDelete = [];
    
    % Delete references in neighbourhood
    % -----------------------------------
    
    for I = 1 : 3
        if TT(iT,I) ~= -1
            
            % Then in the I reference the triangle confines with another
            % triangle
            TT ( TT(iT,I) , TT(iT,I+6) ) = -1;
            TT ( TT(iT,I) , TT(iT,I+6) + 6 ) = -1;
            B ( TT(iT,I+3) , 3:4 ) = [-1 TT(iT,I)];
            
            %Border became a triangulation limit: so:
            BInfo (TT(iT,I+3),2) = 1;
            encroached_borders_add ( TT(iT,I+3) );
            
        else
            % The triangle already is adjacent to outer region. A Border
            % (and may be some vertex) may be deleted
            nBorderToDelete = nBorderToDelete + 1;
            BorderToDelete(nBorderToDelete) =  TT ( iT , I+3 );
            
            % Check if there are also verteces to be deleted
            
            for J = 1:2
                
                jV = B(BorderToDelete(nBorderToDelete),J);
                
                % Find where the border is memorized in VB(jV).B
                for j = 1 : VB(jV).n
                    if VB(jV).B(j) == BorderToDelete(nBorderToDelete)
                        break
                    end
                end
                    
                % Update VB(jV)
                [VB(jV).V] = array_remove_component (VB(jV).V,VB(jV).n,j);
                [VB(jV).B] = array_remove_component (VB(jV).B,VB(jV).n,j);
                VB(jV).n = VB(jV).n - 1;

                if VB ( jV ). n == 0
                    % The vertex must be probably deleted
                    Flag = false;
                    for k = VertexToDelete
                        if jV == k
                            Flag = true;
                            break
                        end
                    end
                    if Flag == false
                        nVertexToDelete = nVertexToDelete + 1;
                        VertexToDelete (nVertexToDelete) = jV;
                    end
              
                end
            end
            
        end
        
    end
    
    % Fill hole in triangles list with TMoved
    % ----------------------------------------
    
    if TMoved == 0
    
        % Triangle is the last in the list
        nT = nT - 1;
        
    else
        
        % Modify TV,TT,TInfo,B
        TV(iT,:) = TV(TMoved,:);
        TT(iT,:) = TT(TMoved,:);
        TInfo(iT) = TInfo(TMoved);
        for I = 1:3
            if TT(TMoved , I) ~= -1
                TT( TT(TMoved,I) , TT(TMoved,I+6) ) = iT;
            end
        end
        B ( TT(iT,4) , 3:4 ) = [ TT(iT,1) , iT ];
        B ( TT(iT,5) , 3:4 ) = [ TT(iT,2) , iT ];
        B ( TT(iT,6) , 3:4 ) = [ TT(iT,3) , iT ];
        
        % Now triangle is the last in the list
        nT = nT - 1;
        
    end
        
    % Delete verteces and fill holes finding VMoved
    % ----------------------------------------------
    
    % Sort VertexToDelete
    if length (VertexToDelete) == 2
        if VertexToDelete(2) > VertexToDelete(1)
            TempSort = VertexToDelete(1);
            VertexToDelete(1) = VertexToDelete(2);
            VertexToDelete(2) = TempSort;
        end
    end

    
    VertexToDeleteCopy = VertexToDelete; % To comparison verteces
    
    for jV = VertexToDelete
        
        if jV == nV
            
            clear TempVertex
            TempVertex = VertexToDeleteCopy (2:length(VertexToDeleteCopy));
            clear VertexToDeleteCopy;
            VertexToDeleteCopy = TempVertex;
            
            %Vertex is the last in the list
            nV = nV - 1;
            
        else
            
            %Find VMoved, vertex that will fill the hole
%%%% CORREZIONE ORTALI             
%            Found == false;
            Found = false;
            VMoved = nV;
            while Found == false
                
                VIsEqualToAnotherV = false;
                
                for j = VertexToDeleteCopy
                    if j == VMoved
                        VMoved = VMoved -1;
                        VIsEqualToAnotherV = true;
                        break
                    end
                end %j
                
                if  VIsEqualToAnotherV == false;
                    Found = true;
                end
                
            end
            
            clear TempVertex
            TempVertex = VertexToDeleteCopy (2:length(VertexToDeleteCopy));
            clear VertexToDeleteCopy;
            VertexToDeleteCopy = TempVertex;
            
            % Modify V
            V (jV,:) = V (VMoved,:);
            VB (jV) = VB (VMoved);
            
            % Modify VB, finding which triangles are to update
            TriangleToUpdate = [];
            nTriangleToUpdate = 0;
            
            for j = 1 : VB (jV).n
                
                OtherVertex = VB (jV).V(j);
                ThisBorder = VB (jV).B(j);
                
                % Add triangles that will be updated
                for I = [3,4];
                
                    Inserted = false;
                    for kT = TriangleToUpdate
                        if B(ThisBorder,I) == kT | B(ThisBorder,I) == -1
                            Inserted = true;
                            break
                        end
                    end
                    
                    if Inserted == false & B(ThisBorder,I) > 0
                        nTriangleToUpdate = nTriangleToUpdate + 1;
                        TriangleToUpdate(nTriangleToUpdate) = B(ThisBorder,I);
                    end
                    
                end
                
                %Update VB and B
                for k = 1: VB(OtherVertex).n
                    if VB(OtherVertex).V(k) == VMoved
                        VB(OtherVertex).V(k) = jV;
                        if B( VB(OtherVertex).B(k) , 1 ) == VMoved
                            B( VB(OtherVertex).B(k) , 1 ) = jV;
                        else
                            B( VB(OtherVertex).B(k) , 2 ) = jV;
                        end
                        break
                    end
                end

                
            end
            
            % Update TV
            for kT = TriangleToUpdate
                for I = [1:3]
                    if TV(kT,I) == VMoved
                        TV(kT,I) = jV;
                        break
                    end
                end
            end
            
            % Update nV
            nV = nV - 1;
            
        end
        
    end %jV
    
    
    % Delete borders and fill holes finding BMoved
    % ----------------------------------------------
    
    % Sort BorderToDelete
    if length (BorderToDelete) == 2
        if BorderToDelete(2) > BorderToDelete(1)
            Temp = BorderToDelete(1);
            BorderToDelete(1) = BorderToDelete(2);
            BorderToDelete(2) = Temp;
        end
    elseif length (BorderToDelete) == 3
        if BorderToDelete(3) > BorderToDelete(1)
            TempSort = BorderToDelete(1);
            BorderToDelete(1) = BorderToDelete(3);
            BorderToDelete(3) = TempSort;
        end
        if BorderToDelete(2) > BorderToDelete(1)
            TempSort = BorderToDelete(1);
            BorderToDelete(1) = BorderToDelete(2);
            BorderToDelete(2) = TempSort;
        end
        if BorderToDelete(3) > BorderToDelete(2)
            TempSort = BorderToDelete(2);
            BorderToDelete(2) = BorderToDelete(3);
            BorderToDelete(3) = TempSort;
        end
    end
    
    BorderToDeleteCopy = BorderToDelete; % To comparison borders
    
    for jB = BorderToDelete
        
        % Update encroachable list
        encroached_borders_remove (jB);
        
        if jB >= nB

            clear TempBorder
            Temp = BorderToDeleteCopy(2:length(BorderToDeleteCopy));
            clear BorderToDeleteCopy;
            BorderToDeleteCopy = Temp;

            %Border is the last in the list
            nB = nB - 1;

            
        else
            
            %Find BMoved, border that will fill the hole
            Found = false;
            BMoved = nB;
            while Found == false
                
                BIsEqualToAnotherB = false;
                
                for j = BorderToDeleteCopy
                    if j == BMoved
                        BMoved = BMoved -1;
                        BIsEqualToAnotherB = true;
                        break
                    end
                end %j
                
                if  BIsEqualToAnotherB == false;
                    Found = true;
                end
                
            end
            
            clear TempBorder
            TempBorder = BorderToDeleteCopy(2:length(BorderToDeleteCopy));
            clear BorderToDeleteCopy;
            BorderToDeleteCopy = TempBorder;
            
            % Find info about borderlist
            BMovedEncroachStatus = BInfo (BMoved,2);
                              
            % Modify B, BInfo
            B (jB,:) = B (BMoved,:);
            BInfo (jB,:) = BInfo (BMoved,:);
            
            % Modify VB
            for I = [1 2]
                
                EndingVertex = B(jB,I);
                
                for k = 1 : VB(EndingVertex).n
                    if VB(EndingVertex).B(k) == BMoved
                        VB(EndingVertex).B(k) = jB;
                        break
                    end
                end
                
            end

            % Modify TT
            for J = 3:4
                if B(BMoved,J) ~= -1
                    for I = 4:6
                        if TT ( B(BMoved,J) , I ) == BMoved
                            TT ( B(BMoved,J) , I ) = jB;
                            break
                        end
                    end
                end
            end
            
            nB = nB - 1;
            
            % Update Encroachable List
            if BMovedEncroachStatus == 1
                encroached_borders_remove (BMoved)
                
                % Modify BCircle (must be downstream with respect to EB update)
                BCircle (jB) = BCircle (BMoved);

                encroached_borders_add (jB)
            else
                
                % Modify BCircle (must be downstream with respect to EB update)
                BCircle (jB) = BCircle (BMoved);
            end

            
        end
               
    end % jB
    
    % delete another triangle
    i = i - 1;
    
end

return
