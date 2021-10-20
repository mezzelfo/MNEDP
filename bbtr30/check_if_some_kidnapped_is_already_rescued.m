function [] = check_if_some_kidnapped_is_already_rescued (Vertex);
%
% [] = try_rescuing_kidnapped_borders ();
%
% This function checks if one of the borders starting form Vertex is a
% subsegment of a KidnappedBorder. The function checks also if the
% remaining part of the segment is another KidnappedBorder or if it
% contains another part of a kidnapped borders
%
% Input of this function is Vertex, the vertex form which the control
% starts


global nKidnappedBorders KidnappedBorders
global V VB BInfo

iKB = 1;
FoundVertex = 0;
OtherVertex = 0;
VertexToControl = [];
while iKB <= nKidnappedBorders
    
    % Find which is OtherVertex (left 0 if not found)
    if KidnappedBorders(iKB,1) == Vertex
        OtherVertex = KidnappedBorders(iKB,2);
    elseif KidnappedBorders(iKB,2) == Vertex
        OtherVertex = KidnappedBorders(iKB,1);
    end
    
    if OtherVertex ~= 0
        
        for iV = 1:VB(Vertex).n;
            
            % Check if Vertex, OtherVertex and ThirdVertex are aligned
            Position = find_alignment_of_memorized_vertex (Vertex,OtherVertex,VB(Vertex).V(iV));
            if Position == 0
                
                ThirdVertex = VB(Vertex).V(iV);
                
                % Check if Vertex & OtherVertex are not aligned in
                % x-direction
                if V(Vertex,1) == V(OtherVertex,1)
                    I = 2;
                else
                    I = 1;
                end
                
                % Check if ThirdVertex is between Vertex & OtherVertex
                if sign(V(Vertex,I)-V(ThirdVertex,I)) == sign(V(ThirdVertex,I)-V(OtherVertex,I))
                    
                    % A KidnappedBorder is found!
                    
                    % Modify BorderInfo and KidnappingBorder
                    BorderInfo(VB(Vertex).B(iV),1) = KidnappedBorders(iKB,4);
                    BorderInfo(VB(Vertex).B(iV),3) = KidnappedBorders(iKB,5);
                    KidnappedBorders(iKB,1:2) = [ThirdVertex OtherVertex];
                                        
                    % Check if the newformed border already exists in current triangulation
                    Found = false;
                    for j = 1 : VB(ThirdVertex).n
                        if VB(ThirdVertex).V(j) == OtherVertex
                            % Another kidnapped border (and already rescued) is found
                            Found = true;
                            BorderInfo(VB(ThirdVertex).B(j),1) = KidnappedBorders(iKB,4);
                            BorderInfo(VB(ThirdVertex).B(j),3) = KidnappedBorders(iKB,5);
                            
                            %%% KidnappedBorder(1:iKB-1) =
                            %%% KidnappedBorder(1:iKB-1); that's implicit
                            KidnappedBorder (iKB:nKidnappedBorders-1) = KidnappedBorder (iKB+1:nKidnappedBorders);
                            nKidnappedBorders = nKidnappedBorders - 1;
                            iKB = iKB - 1; %Because KidnBord(iKB) has been deleted, and other controls must be made
                            break
                        end
                    end
                                        
                    % Check if other KidnappedBorder may be found,
                    if Found == false
                        
                        FoundVertex = FoundVertex + 1;
                        VertexToControl (FoundVertex) = ThirdVertex;
                        
                    end    
                end
                
            end
            
        end
        
    end
        
    iKB = iKB + 1;
    
    % Recursively controls newformed domain
    for CurrentVertex = VertexToControl
        check_if_some_kidnapped_is_already_rescued (CurrentVertex);
    end
     
end

return