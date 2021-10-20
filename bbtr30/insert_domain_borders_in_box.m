function [] = insert_domain_borders_in_box ();
%
% [] = insert_domain_borders_in_box ();
%
% This function inserts every kidnapped border (* see attachments) in the
% triangulation created in the CB.

global nKidnappedBorders

while nKidnappedBorders ~= 0

    try_rescuing_kidnapped_borders
   
end

return