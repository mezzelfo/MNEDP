function [Length] = find_length_of_memorized_border (Border);
%
% [Length] = find_length_of_memorized_border (Border);
%
% This function find the length of a memorized border
%
% Input & Output of this function are:
%
% Border: the border whoselength is unknown
% Length: the length of Border
%

global B V

Length = sqrt( ( V(B(Border,1),1) - V(B(Border,2),1) )^2 + ( V(B(Border,1),2) - V(B(Border,2),2) )^2 );

return