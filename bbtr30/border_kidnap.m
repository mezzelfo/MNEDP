function [] = border_kidnap (Border,Triangle);
%
% [] = border_kidnap (Border);
%
% This function kidnaps Border (*see attachments), updating necessary list
%
% Input of this function are:
%
% Border: The Border to kidnap
% Triangle: a triangle near to Border (may be useful later to rescue it) 
%

 
global B BInfo
global nKidnappedBorders KidnappedBorders

nKidnappedBorders = nKidnappedBorders + 1;
KidnappedBorders (nKidnappedBorders,:) = [B(Border,1:2) Triangle BInfo(Border,1) BInfo(Border,3)]; 

return