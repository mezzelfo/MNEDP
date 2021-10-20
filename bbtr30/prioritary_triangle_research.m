function [T,RelativePosition] = prioritary_triangle_research (xV,yV,Tstart);
%
% [T RelativePosition] = prioritary_triangle_research (xV,yV,Tstart)
%
% This function finds which triangle enclosures vertex (xV,yV). The
% research starts from triangle Tstart, and then if the vertex doesn't lyes
% inside it the focus is moved to his neighbour triangles. The procedure (PTR)
% is applied again until the triangle is found or until the focused triangle
% too "far" from Tstart (* see attachments to find other informations)
%
% Inputs & Output of this function are
%
% T: the triangle that enclosure Vertex (xV,yV)
% RelativePosition: a marker that says if (xV,yV) lyes on a border (gives the
%                   reference to that border) or inside T (gives back value
%                   0);
%
% (xV,yV): the vertex to be found in the grid;
% Tstart: the triangle from which the reserach starts.
%

global GenerationUntilReserach

% Check if vertex lyes in Tstart
RelativePosition = check_if_in_triangle (xV,yV,Tstart);

if RelativePosition >= 0
    
    T = Tstart;
    
else
    
    % Check if vertex lyes in TstartNeighbours
    global TT
    
    for CurrentTriangle = TT(Tstart,1:3)
        
        if CurrentTriangle ~= -1
            RelativePosition = check_if_in_triangle (xV,yV,CurrentTriangle);
            if RelativePosition >= 0
                T = CurrentTriangle;
                break
            end
        end
        
    end

    % Reserch with PTR

    if RelativePosition < 0
        
        global nT
        global GenerationUntilResearch
        
        % Update CheckedTriangle list
        CheckedTriangle = zeros (nT,1); % 0 --> Unchecked; 1 --> in check list; 2--> Checked
        
        CheckedTriangle (Tstart) = 2;
        if TT(Tstart,1) ~= -1
            CheckedTriangle (TT(Tstart,1)) = 2;
        end
        if TT(Tstart,2) ~= -1
            CheckedTriangle (TT(Tstart,2)) = 2;
        end
        if TT(Tstart,3) ~= -1
            CheckedTriangle (TT(Tstart,3)) = 2;
        end

        % Create first TriangleToCheck list
        nOldTriangleToCheck = 0;
        nNewTriangleToCheck = 0;
        OldTriangleToCheck = zeros (min(3^GenerationUntilResearch,nT),1); %It's a conservative extimation
        NewTriangleToCheck = zeros (min(3^GenerationUntilResearch,nT),1); %It's a conservative extimation

        for CurrentTriangle = TT(Tstart,1:3)
            for I = 1:3
                if CurrentTriangle ~= -1
                    if TT(CurrentTriangle,I) ~= -1
                        if CheckedTriangle ( TT(CurrentTriangle,I) ) == 0
                            CheckedTriangle ( TT(CurrentTriangle,I) ) = 1; %In checklist
                            nOldTriangleToCheck = nOldTriangleToCheck + 1;
                            OldTriangleToCheck (nOldTriangleToCheck) = TT(CurrentTriangle,I);
                        end
                    end
                end
            end
        end
        
        % Launch Real PTR
        CurrentGeneration = 1;
        
        while (RelativePosition < 0) & (CurrentGeneration <= GenerationUntilResearch)
            
            CurrentGeneration = CurrentGeneration + 1;
            
            for iCT = 1:nOldTriangleToCheck
                
                CurrentTriangle = OldTriangleToCheck (iCT);
                % Checks if vertex lyes in CurrentTriangle
                RelativePosition = check_if_in_triangle (xV,yV,CurrentTriangle);
                CheckedTriangle (CurrentTriangle) = 2;
                
                if RelativePosition >= 0
                    % Vertex lyes inside triangle!
                    T = CurrentTriangle;
                    break
                else
                    % Vertex Lyes outside triangle.
                    for I = 1:3
                        if TT(CurrentTriangle,I) ~= -1
                            if CheckedTriangle ( TT(CurrentTriangle,I) ) == 0
                                CheckedTriangle ( TT(CurrentTriangle,I) ) = 1; %In checklist
                                nNewTriangleToCheck = nNewTriangleToCheck + 1;
                                NewTriangleToCheck (nNewTriangleToCheck) = TT(CurrentTriangle,I);
                            end
                        end
                    end
                end
                
            end
            
            nOldTriangleToCheck = nNewTriangleToCheck;
            OldTriangleToCheck(1:nOldTriangleToCheck) = NewTriangleToCheck (1:nOldTriangleToCheck);
            nNewTriangleToCheck = 0;
            
        end
        
        if RelativePosition < 0
            
            % Unsorted & unoptimized reserarch in the whole grid
            
            for CurrentTriangle = 1 : nT
                
                if CheckedTriangle (CurrentTriangle) == 0
                    
                    RelativePosition = check_if_in_triangle (xV,yV,CurrentTriangle);
                    
                    if RelativePosition >= 0
                        % Vertex lyes inside triangle!
                        T = CurrentTriangle;
                        break
                    end
                    
                end
                
            end

        end

        
    end

    
end

if RelativePosition < 0
    disp ('WARNING: Point outside the domain is found form prioritary_triangle_research funcion')
end

return
            





