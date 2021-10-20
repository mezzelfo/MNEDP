function [SortedInsertionList,NextMaxRoof,WhatHasBeenRefined] = create_inserction_list (RefiningOptions,RefiningMaxRoof,PreviousMaxRoof,WhatHasBeenRefined,WhereWeMustArriveInRefining);
%
% [SortedInsertionList,NextMaxRoof,WhatHasBeenRefined] = create_inserction_list (RefiningOptions,PreviousMaxRoof,PreviousRoof,WhatHasBeenRefined,WhereWeMustArriveInRefining);
%
% This function creates OrderedInserctionList, a matrix containing the list
% of next triangles to refine and determines what must be refined.
%
% Ouput & Input of this function are:
%
% RefiningOptions: a structure that explain how to impose the
%   condition that will identify the grid (i.e. the minumum angle checked
%   and the maximum area allowed). It is possible to refine only a
%   particular subregion of the domain. (Particular usages of this variable
%   may be found in the attachments
% PreviousMaxRoof: says the theorical after which refining_grid should
%   has stopped insercting triangles
% PreviousRoof: says after how many triangles previously refining_grid
%   has stopped inserction
% WhatHasBeenRefined: is a variable that indicates what the program is
%   refining. It may assums following values:
%   -1: working on area of domain - still some elements to split
%   1: working on area of domain - not any more element to split
%   -2: working on quality of domain - still some elements to split
%   2: working on quality of domain - not any more element to split
%   -(2i+1), with i>1: working on area of i refining subregion - still some elements to split
%   (2i+1), with i>1: working on area of i refining subregion - not any more element to split
%   -(2i+2), with i>1: working on quality of i refining subregion - still some elements to split
%   (2i+2), with i>1: working on quality of i refining subregion - not any more element to split
% WhereWeMustArriveInRefining: is a variable that indicates where the
%   program have to stop in refining. It may assums following values
%   0: no refining required.
%   1: last refining is on area of domain.
%   2: last refining is on quality of domain.
%   (2i+1), with i>1: last refining is on area of i refining subregion.
%   (2i+2), with i>1: last refining is on quality of i refining subregion.
%
% SortedInserctionList: an ordered matrix containing in the first column a
%   reference to which triangles must be inserted, and in the second his
%   area. This is useful because is the only way to know if later will not
%   be refined an accettable triangle.
% PreviousMaxRoof: says the maximum threshold after which refining_grid will
%   stop insercting triangles
%

global FirstMaxRefiningRoof IncreasingRefiningRoof RefiningFolderName
global TInfo
global nT

% First of all find how will be next roof and what will be refined
% -----------------------------------------------------------------

if WhatHasBeenRefined == 0 
    
    % This is the first refining. The value must be chosen with experience aid
    NextMaxRoof = min( ceil(FirstMaxRefiningRoof * nT),10);
    
    % Determines WhatHasBeenRefined for the first time
    if RefiningOptions.CheckArea == 'y' | RefiningOptions.CheckArea == 'Y'
        WhatHasBeenRefined = -1;
    elseif RefiningOptions.CheckAngle == 'y' | RefiningOptions.CheckAngle == 'Y'
        WhatHasBeenRefined = -2;
    elseif RefiningOptions.Subregions.Subregion(1).CheckArea == 'y' | RefiningOptions.Subregions.Subregion(1).CheckArea == 'Y'
        WhatHasBeenRefined = -3;
    else
        WhatHasBeenRefined = -4;
    end
    
elseif WhatHasBeenRefined > 0
    
    % An inserction kind (i.e.: area in main region, angle in subregion 2 )
    % has ende. The value must be chosen with experience aid
    NextMaxRoof = min( ceil(FirstMaxRefiningRoof * nT),10);
    
    % Updates WhatHasBeenRefined
    if WhatHasBeenRefined == 1
        
        % Check angle of main region
        if RefiningOptions.CheckAngle == 'y' | RefiningOptions.CheckAngle == 'Y'
            WhatHasBeenRefined = -2;
        elseif RefiningOptions.Subregions.Subregion(1).CheckArea == 'y' | RefiningOptions.Subregions.Subregion(1).CheckArea == 'Y'
            WhatHasBeenRefined = -3;
        else
            WhatHasBeenRefined = -4;
        end
        
    else
        
        % A generalized control may be done
        if mod (WhatHasBeenRefined,2) == 0
            
            % We finished an Angle control
            NextRegion = (WhatHasBeenRefined - 2)/2 + 1;
            if RefiningOptions.Subregions.Subregion(NextRegion).CheckArea == 'y' | RefiningOptions.Subregions.Subregion(NextRegion).CheckArea == 'Y'
                WhatHasBeenRefined = - (WhatHasBeenRefined + 1);
            elseif RefiningOptions.Subregions.Subregion(NextRegion).CheckAngle == 'y' | RefiningOptions.Subregions.Subregion(NextRegion).CheckAngle == 'Y'
                WhatHasBeenRefined = - (WhatHasBeenRefined + 2);
            else
                WhatHasBeenRefined = - (WhatHasBeenRefined + 3);
            end
            
        else
            
            % We finished an Area control
            CurrentRegion = (WhatHasBeenRefined - 1)/2;
            if RefiningOptions.Subregions.Subregion(CurrentRegion).CheckAngle == 'y' | RefiningOptions.Subregions.Subregion(CurrentRegion).CheckAngle == 'Y'
                WhatHasBeenRefined = - (WhatHasBeenRefined + 1);
            elseif RefiningOptions.Subregions.Subregion(CurrentRegion + 1).CheckArea == 'y' | RefiningOptions.Subregions.Subregion(CurrentRegion +1).CheckArea == 'Y'
                WhatHasBeenRefined = - (WhatHasBeenRefined + 2);
            else
                WhatHasBeenRefined = - (WhatHasBeenRefined + 3);
            end
            
        end
        
    end

elseif  WhatHasBeenRefined < 0
    
    % WhatHasBeenRefined is manteined. Only NextMaxRoof is modified.
    if PreviousMaxRoof == RefiningMaxRoof
        
        NextMaxRoof = IncreasingRefiningRoof * PreviousMaxRoof;
        
    else
        
        NextMaxRoof = PreviousMaxRoof;
        
    end
        
end

% Create InserctionList
% ----------------------

iT = zeros (nT,2); %IMP: iT HERE IS A MATRIX, NOT A SCALAR!
niT = 0;
NumberOfTriangles = 0;
LimitReached = true; %This variable remains true if there are no more elements to split

if mod (WhatHasBeenRefined,2) == 1
    
    % This is an area control
    
    if WhatHasBeenRefined == -1
        
        AreaThreshold = RefiningOptions.AreaValue;
        
        for jT = 1: nT
            
            if TInfo(jT).Area > AreaThreshold

                NumberOfTriangles = NumberOfTriangles + 1;
                
                if NumberOfTriangles > NextMaxRoof
                    LimitReached = false;
                    break
                else
                    niT = niT + 1;
                    iT(niT,:) = [jT TInfo(jT).Area];
                end
                
                
            end
            
        end
        
    else
        
        eval (['cd ',RefiningFolderName]);
        Region = -(WhatHasBeenRefined + 1)/2;
        AreaThreshold = RefiningOptions.Subregions.Subregion(Region).AreaValue;
        
        for jT = 1: nT
            
            if TInfo(jT).Area > AreaThreshold
                
                % Checks also if jT is in region
                TinRegion = feval ( RefiningOptions.Subregions.Subregion(Region).Type , RefiningOptions.Subregions.Subregion(Region).Descriptors, TInfo(jT).CG );
                if TinRegion == true
                    
                    NumberOfTriangles = NumberOfTriangles + 1;
                    
                    if NumberOfTriangles > NextMaxRoof
                        LimitReached = false;
                        break
                    else
                        niT = niT + 1;
                        iT(niT,:) = [jT TInfo(jT).Area];
                    end

                end
                
            end
            
        end
        
        cd ..
        
    end
    
else

    % This is a quality control
    
    if WhatHasBeenRefined == -2
        
        BThreshold = 1 / (2 * sin (RefiningOptions.AngleValue*pi/180) );
        
        for jT = 1: nT
            
            if TInfo(jT).B > BThreshold
                
                NumberOfTriangles = NumberOfTriangles + 1;
                
                if NumberOfTriangles > NextMaxRoof
                    LimitReached = false;
                    break
                else
                    niT = niT + 1;
                    iT(niT,:) = [jT TInfo(jT).B];
                end
                
            end
            
        end
        
    else
        
        eval (['cd ',RefiningFolderName]);
        Region = -(WhatHasBeenRefined + 2)/2;
        BThreshold = 1 / (2 * sin (RefiningOptions.Subregions.Subregion(Region).AngleValue*pi/180) );
        
        for jT = 1: nT
            
            if TInfo(jT).B > BThreshold
                
                % Checks also if jT is in region
                TinRegion = feval ( RefiningOptions.Subregions.Subregion(Region).Type , RefiningOptions.Subregions.Subregion(Region).Descriptors, TInfo(jT).CG );
                if TinRegion == true
                    
                    NumberOfTriangles = NumberOfTriangles + 1;
                
                    if NumberOfTriangles > NextMaxRoof
                        LimitReached = false;
                        break
                    else
                        niT = niT + 1;
                        iT(niT,:) = [jT TInfo(jT).B];
                    end
                    
                end
                
            end
            
        end
        cd ..;

    end
    
end

% Create SortedInserctionList, sorting iT
% ----------------------------------------

Temp = iT(1:niT,:);
Temp = my_sort_descending (Temp,niT,2,2);

if NextMaxRoof <= niT
    SortedInsertionList = Temp(1:NextMaxRoof,:);
else
    SortedInsertionList = Temp;
    NextMaxRoof = niT;
end


if LimitReached == true
    WhatHasBeenRefined = abs(WhatHasBeenRefined);
end
    
return