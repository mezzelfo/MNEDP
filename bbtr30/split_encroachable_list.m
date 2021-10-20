function [] = split_encroachable_list (CurrentDomain);
%
% [] = split_encroachable_list (CurrentDomain);
%
% This function splits a subdomain assigning every border to a specific
% splitted subdomain. The split is horizontal (index = 1) or vertical
% (index = 2) according to the larger dimension of domain to be splitted
% 
% Input of this function is:
%
% CurrentDomain : the domain to be splitted
%

global EB
global BCircle

% Set two next subdomains

SubDomainLower = EB.nBr + 1;
SubDomainUpper = EB.nBr + 2;
EB.nBr = EB.nBr + 2;
EB.Br(CurrentDomain).Unsplitted = false;

EB.Tree (CurrentDomain,1:2) = [SubDomainLower SubDomainUpper];
EB.Tree (SubDomainLower:SubDomainUpper,1:2) = [0 0; 0 0];

for iD = [SubDomainLower SubDomainUpper]
    EB.Br(iD).Unsplitted = true;
    EB.Br(iD).Leaves = [];
    EB.Br(iD).nLeaves = 0;
end



% Determines if next splitting will be horizontal or vertical

DomainLimits = EB.Br(CurrentDomain).SplittingLimits;
if abs (DomainLimits(1) - DomainLimits(2)) > abs (DomainLimits(3) - DomainLimits(4))
    
    EB.Br(CurrentDomain).KindOfSplitting = 1;
    EB.Br(CurrentDomain).SplittingThreshold = (DomainLimits(1) + DomainLimits(2))/2;
    
    LowerDomainLimits = [EB.Br(CurrentDomain).SplittingThreshold,DomainLimits(2),DomainLimits(3),DomainLimits(4) ];
    UpperDomainLimits = [DomainLimits(1),EB.Br(CurrentDomain).SplittingThreshold,DomainLimits(3),DomainLimits(4) ];
    
    EB.Br(SubDomainLower).SplittingLimits = LowerDomainLimits;
    EB.Br(SubDomainUpper).SplittingLimits = UpperDomainLimits;

else
    
    EB.Br(CurrentDomain).KindOfSplitting = 2;
    EB.Br(CurrentDomain).SplittingThreshold = (DomainLimits(3) + DomainLimits(4))/2;
    
    LowerDomainLimits = [DomainLimits(1),DomainLimits(2),EB.Br(CurrentDomain).SplittingThreshold,DomainLimits(4)];
    UpperDomainLimits = [DomainLimits(1),DomainLimits(2),DomainLimits(3),EB.Br(CurrentDomain).SplittingThreshold];

    EB.Br(SubDomainLower).SplittingLimits = LowerDomainLimits;
    EB.Br(SubDomainUpper).SplittingLimits = UpperDomainLimits;

end



% Assign every border to a specific subdomain (or to both)

for iB = EB.Br(CurrentDomain).Leaves
    
    % Find circumcenter and circumradius
    xC = BCircle(iB).Circumcenter(1);
    yC = BCircle(iB).Circumcenter(2);
    r = sqrt ( BCircle(iB).r2 );
    
    % Find if iB belongs to UpperDomainLimit
    if ( xC+r >= UpperDomainLimits(2)) & ( xC-r <= UpperDomainLimits(1)) & ( yC+r >= UpperDomainLimits(4)) & ( yC-r <= UpperDomainLimits(3))
        
        EB.Br(SubDomainUpper).nLeaves = EB.Br(SubDomainUpper).nLeaves + 1;
        EB.Br(SubDomainUpper).Leaves ( EB.Br(SubDomainUpper).nLeaves ) = iB;
        
    end
    
    % Find if iB belongs to LowerDomainLimit
    if ( xC+r >= LowerDomainLimits(2)) & ( xC-r <= LowerDomainLimits(1)) & ( yC+r >= LowerDomainLimits(4)) & ( yC-r <= LowerDomainLimits(3))
        
        EB.Br(SubDomainLower).nLeaves = EB.Br(SubDomainLower).nLeaves + 1;
        EB.Br(SubDomainLower).Leaves ( EB.Br(SubDomainLower).nLeaves ) = iB;
        
    end
    
end

% Delete references in splitted domain
EB.Br(CurrentDomain).nLeaves = 0;
EB.Br(CurrentDomain).Leaves = [];

% Eventually splits newborn domain
if EB.Br(SubDomainLower).nLeaves == 0
    split_encroachable_list (SubDomainUpper);
elseif EB.Br(SubDomainUpper).nLeaves == 0
    split_encroachable_list (SubDomainLower);
end
    
    
return    