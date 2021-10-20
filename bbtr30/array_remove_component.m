function [OutputArray] = array_remove_component (InputArray,n,i);
%
% [OutputArray] = array_remove_component (InputArray,n,i)
%
% This function erase the component in i-place from the array InputArray.
%
% Output & Input of this function are:
%
% InputArray: the array to be shortened
% n: the size of InputArray
% i: the position of the element that must be removed
%
% OutputArray: the shortened array
%

OutputArray(1:i-1) = InputArray(1:i-1);
OutputArray(i:n-1) = InputArray(i+1:n);

return