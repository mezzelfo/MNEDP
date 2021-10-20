function [SubDomain,nSubDomain] = encroached_borders_search_in_depth (CurrentDomain,Border,SubDomain,nSubDomain);
%
% [SubDomain,nSubDomain] = encroached_borders_search_in_depth (CurrentDomain,Border,SubDomain,nSubDomain);
%
% This functions discovers to which existing subsets of domains belongs
% Border and puts the risults in SubDomain
%
% Input & Outputs of this function are:
%
% CurrentDomain : the domain from which start the research
% Border: the Border to be searched
% SubDomain: the list of subdomains to which we already kows that Border
%   belongs to. It is updated by the function
% nSubDomain: the length of SubDomain
%

global EB
global BCircle

Found = false;
while Found == false
    
    if EB.Br ( CurrentDomain ).Unsplitted == true
        
        % A subdomain is found!
        nSubDomain = nSubDomain + 1;
        SubDomain(nSubDomain) = CurrentDomain;
        Found = true;
        
    else
        
        % Search in subdomains
        I = EB.Br ( CurrentDomain ).KindOfSplitting;
        
        if ( BCircle(Border).Circumcenter(I) - sqrt (BCircle(Border).r2) ) <= EB.Br ( CurrentDomain ).SplittingThreshold;
            
            if ( BCircle(Border).Circumcenter(I) + sqrt (BCircle(Border).r2) ) >= EB.Br ( CurrentDomain ).SplittingThreshold;
                
                % A double subdomain is found
                [SubDomain,nSubDomain] = encroached_borders_search_in_depth ( EB.Tree(CurrentDomain,1) ,Border,SubDomain,nSubDomain);
                NextDomain = EB.Tree(CurrentDomain,2);
                
            else
                
                NextDomain = EB.Tree(CurrentDomain,1);
                
            end
            
        else
            
            NextDomain = EB.Tree(CurrentDomain,2);
            
        end
        
        CurrentDomain = NextDomain;
        
    end
    
end

return
        