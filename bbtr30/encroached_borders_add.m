function [] = encroached_borders_add(Border);
%
% [] = encroached_borders_add(Border);
%
% This function adds Border to Encroachable list (*see attachments for
% further informations).
%
% Input of this funcion is Border, the border to add to the list

global EB
global SplittingBordersLimits

% Search to what subdomains the border belongs to
nSubDomain = 0;
SubDomain = [];
CurrentDomain = 1;

[SubDomain nSubDomain] = encroached_borders_search_in_depth (CurrentDomain,Border,SubDomain,nSubDomain);

for i = 1:nSubDomain
    
    iS = SubDomain (i);
    EB.Br( iS ).nLeaves = EB.Br( iS ).nLeaves + 1;
    EB.Br( iS ).Leaves ( EB.Br( iS ).nLeaves ) = Border; %EB.Tree(CurrentSubDomain,3)
    
    if EB.Br( iS ).nLeaves >= SplittingBordersLimits
        
        split_encroachable_list (iS);
        
    end
    
end

return