function [nTApprox,nVApprox,nBApprox] = first_size_variables_extimate (nInputVertex);
%
% [nTApprox,nVApprox,nBApprox] = first_size_variables_extimate (nInputVertex);
%
% This function gives a first extimate of the size of the 
% elements of the grid 
%
% Inputs & Outputs of this function are:
%
% nInputVertex : the number of input verteces
% nTApprox :  the extimated number of triangle necessary to insert domain
% nVApprox :  the extimated number of verteces necessary to insert domain
% nBApprox :  the extimated number of borders necessary to insert domain
%

global FirstExtimateTmultiplier FirstExtimateBmultiplier FirstExtimateVmultiplier

nTApprox = nInputVertex * FirstExtimateTmultiplier;
nVApprox = nInputVertex * FirstExtimateVmultiplier;
nBApprox = nInputVertex * FirstExtimateBmultiplier;

return