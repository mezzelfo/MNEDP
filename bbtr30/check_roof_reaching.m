function WhatHasBeenRefined = check_roof_reaching (WhatHasBeenRefined,RefiningOptions);
%
% [WhatHasBeenRefined] = check_roof_reaching (WhatHasBeenRefined,RefiningOptions);
%
% This function checks if every triangle of the grid satisfy the condition
% that WhatHasBeenRefined suggests that is finished. If this event doesn't
% occur the function notes it the output
%
% Output & input of this function are:
%
% RefiningOptions: a structure that explain how to impose the
%   condition that will identify the grid (i.e. the minumum angle checked
%   and the maximum area allowed). It is possible to refine only a
%   particular subregion of the domain. (Particular usages of this variable
%   may be found in the attachments
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
%

global nT TInfo
global RefiningFolderName;


if WhatHasBeenRefined < 3
    
    %Working on the domain
    
    if WhatHasBeenRefined == 1
        
        % Working on area
        RoofReached = true;
        for iT = 1:nT
            if TInfo(iT).Area > RefiningOptions.AreaValue
                RoofReached = false;
                break
            end
        end
        
    else
        
        % Working on quality
        RoofReached = true;
        BThreshold = 1 / (2 * sin (RefiningOptions.AngleValue*pi/180) );
        for iT = 1:nT
            if TInfo(iT).B > BThreshold
                RoofReached = false;
                break
            end
        end
        
    end
    
else
    
    if mod(WhatHasBeenRefined,2) == 1
        
        % Working on area
        Region = RefiningOptions.Subregions.Subregion( (WhatHasBeenRefined-1)/2 );
        RoofReached = true;
        eval (['cd ',RefiningFolderName]);
        for iT = 1:nT
            if TInfo(iT).Area > Region.AreaValue
                
                % Apply current method
                TinRegion = feval ( Region.Type , Region.Descriptors, TInfo(iT).CG );
                if TinRegion == true
                    RoofReached = false;
                    break
                end
                
            end
        end
        cd ..
        
    else
        
        % Working on quality
        Region = RefiningOptions.Subregions.Subregion( (WhatHasBeenRefined-2)/2 );
        BThreshold = 1 / (2 * sin (Region.AngleValue*pi/180) );
        RoofReached = true;
        eval (['cd ',RefiningFolderName]);
        for iT = 1:nT
            if TInfo(iT).B > BThreshold
                
                % Apply current method
                TinRegion = feval ( Region.Type , Region.Descriptors, TInfo(iT).CG );
                if TinRegion == true
                    RoofReached = false;
                    break
                end
                
            end
        end
        cd ..
        
    end

end

% Checks if the roof is reached

if RoofReached == false
    WhatHasBeenRefined = - WhatHasBeenRefined;
end

return

    