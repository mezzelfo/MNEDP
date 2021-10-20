function [] = check_if_kidnapped_are_found (BorderToExamine);
%
% [] = check_if_kidnapped_are_found (BorderToExamine);
%
% This function checks if BorderToExamine is one of the kidnapped borders (*see
% attachments). If it's true the function rescue (*see
% attachments) the border.
%
% Input of this function is the BorderToExamine that may be a kidnapped one.

global nKidnappedBorders 

for i = 1 : nKidnappedBorders 

    global B BInfo
    global KidnappedBorders

    if B(BorderToExamine,1) == KidnappedBorders(i,1)
        if B(BorderToExamine,2) == KidnappedBorders(i,2)
            
            %Found a kidnapped border! Updating references
            BInfo(BorderToExamine,:) = [KidnappedBorders(i,4) 0 KidnappedBorders(i,5)];
            KidnappedBorders(i:nKidnappedBorders-1 , :) = KidnappedBorders(i+1:nKidnappedBorders , :);
            nKidnappedBorders = nKidnappedBorders - 1;
            break
            
        end    
    elseif B(BorderToExamine,2) == KidnappedBorders(i,1)
        if B(BorderToExamine,1) == KidnappedBorders(i,2)

            %Found a kidnapped border! Updating references
            BInfo(BorderToExamine,:) = [KidnappedBorders(i,4) 0 KidnappedBorders(i,5)];
            KidnappedBorders(i:nKidnappedBorders-1 , :) = KidnappedBorders(i+1:nKidnappedBorders , :);
            nKidnappedBorders = nKidnappedBorders - 1;
            break
            
        end    
    end
end 

return