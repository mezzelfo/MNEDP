function [TriangleKind] = autoidentify_grid();
%
% [TriangleKind] = autoidentify_grid()
%
% This function considers an existing grid containing triangles outside the
% region and, using autoidentify alghoritm (*see attachments) discovers if
% every triangle lyes inside or aoutside the domain.
%
% Output of this function are:
%
% TriangleKind: an array that assumes a value for every triangle in the
%   grid. 1 means that the triangle must be keed, 0 means that the triangle
%   must be deleted
%

global nT
global TT BInfo

% Define TriangleKind
TriangleKind = ones (nT,1); 

%    % TriangleKind will assume these values: 
%    % 1 if triangle is still not assigned to a region
%    % 0 if sure that triangle is outside domain
%    % 2 if sure that triangle is inside the domain
%    % 3 if triangle is studied in current region


% Start with auto_identify_alghoritm    
nMarkedTriangle = 0;    
LastStudiedTriangle = 1;

while nMarkedTriangle < nT
    
    % Define a new region and region variables
    CurrentRegion = [];
    nCurrentRegion = 0;
    InnerFlag = true;
    HoleFlag = true;
    
    % Determine which is first triangle still not studied
    while TriangleKind (LastStudiedTriangle) ~= 1
        LastStudiedTriangle = LastStudiedTriangle + 1;
    end    
    
    % Determine which triangles belong to region
    % --------------------------------------------
    
    nOldAdvancingFront = 1;
    OldAdvancingFront = LastStudiedTriangle;
    
    while nOldAdvancingFront ~= 0
        
        nNewAdvancingFront = 0;
        NewAdvancingFront = [];
        
        for i = 1:nOldAdvancingFront
            iT = OldAdvancingFront (i);
            
            if TriangleKind (iT) == 1
                
                % Mark that iT belongs to CurrentRegion
                TriangleKind(iT) = 3;
                nCurrentRegion = nCurrentRegion + 1;
                CurrentRegion(nCurrentRegion) = iT;
                
                % Study neighbourhood & Trianglekind
                for I = 1:3
                    
                    if TT(iT,I) ~= -1
                        
                        if BInfo( TT(iT,I+3) , 1 ) == 1 | BInfo( TT(iT,I+3) , 1 ) == 3 %(Domain limit or segment limit)
                            
                            HoleFlag = false;
                            
                        elseif TriangleKind (TT(iT,I)) == 1 & BInfo( TT(iT,I+3) , 1 ) == 0 & BInfo( TT(iT,I+3) , 2 ) == 0
                            
                            % Check also next triangle
                            nNewAdvancingFront = nNewAdvancingFront + 1;
                            NewAdvancingFront (nNewAdvancingFront) = TT(iT,I);
                            
                        end
                        
                    else
                        
                        HoleFlag = false;
                        
                        if BInfo( TT(iT,I+3) , 1 ) == 0
                            InnerFlag = false;
                        end   
                        
                    end
                    
                end
                
            end
            
        end
        
        nOldAdvancingFront = nNewAdvancingFront;
        OldAdvancingFront (1:nOldAdvancingFront) = NewAdvancingFront (1:nNewAdvancingFront);

    end


    
    nMarkedTriangle = nMarkedTriangle + nCurrentRegion;
    % --------------------------------------------------
    
    
    % Determine the kind of region
    if ~(InnerFlag == true & HoleFlag == false)
        
        % The region is an outer region
        KindOfRegion = 0;
        
    else
        
        % The region is an inner region
        KindOfRegion = 2;
        
    end

    % Assign every triangle to the region
    for i = CurrentRegion
        TriangleKind (i) = KindOfRegion;
    end
    
end

% Assign to TriangleKind correct meaning:
for i = 1:nT
    
    if TriangleKind(i) == 2
        TriangleKind(i) = 1;
    end
    
end

TriangleKind;

return