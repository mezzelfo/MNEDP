function [] = encroached_borders_remove(Border);
%
% [] = encroached_borders_remove(Border);
%
% This function removes Border to Encroachable list (*see attachments for
% further informations).
%
% Input of this funcion is Border, the border to add to the list

global EB

% Search to what subdomains the border belongs to
nSubDomain = 0;
SubDomain = [];
CurrentDomain = 1;

[SubDomain nSubDomain] = encroached_borders_search_in_depth (CurrentDomain,Border,SubDomain,nSubDomain);

for i = 1:nSubDomain
    
    iS = SubDomain (i);
    
    for iB = 1:EB.Br(iS).nLeaves
        if EB.Br(iS).Leaves(iB) == Border
            break
        end
    end
    
    EB.Br(iS).Leaves = array_remove_component (EB.Br(iS).Leaves,EB.Br(iS).nLeaves,iB);
    EB.Br(iS).nLeaves = EB.Br(iS).nLeaves - 1;
    
end

return    