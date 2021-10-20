function [WhereWeMustArriveInRefining] = find_end_of_refining (RefiningOptions)
%
% [WhereWeMustArriveInRefining] = find_end_of_refining (RefiningOptions);
%
% This function discovers gives a value to the variable
% WhereWeMustArriveInRefining, that is used to discover where stop the
% refining
%
% Output & input of this function are:
%
% RefiningOptions: a structure that explain how to impose the
%   condition that will identify the grid (i.e. the minumum angle checked
%   and the maximum area allowed). It is possible to refine only a
%   particular subregion of the domain. (Particular usages of this variable
%   may be found in the attachments
% WhereWeMustArriveInRefining: is a variable that indicates where the
%   program have to stop in refining. It may assums following values
%   0: no refining required.
%   1: last refining is on area of domain.
%   2: last refining is on quality of domain.
%   (2i+1), with i>1: last refining is on area of i refining subregion.
%   (2i+2), with i>1: last refining is on quality of i refining subregion.
%

if length(RefiningOptions.Subregions) == 0
    
    %There are no refining subregions
    if RefiningOptions.CheckAngle == 'y' | RefiningOptions.CheckAngle == 'Y'
        WhereWeMustArriveInRefining = 2;
    elseif RefiningOptions.CheckArea == 'y' | RefiningOptions.CheckArea == 'Y'
        WhereWeMustArriveInRefining = 1;
    else
        WhereWeMustArriveInRefining = 0;
    end
    
else
    
    % At least a region to refine exists
    i = length (RefiningOptions.Subregions.Subregion);
    if RefiningOptions.Subregions.Subregion(i).CheckAngle == 'y' | RefiningOptions.Subregions.Subregion(i).CheckAngle == 'Y'
        WhereWeMustArriveInRefining = 2*i + 2;
    else
        WhereWeMustArriveInRefining = 2*i + 1;
    end
        
end

return