function [] = reinitialize_variables (newT,newV,newB)
%
% [] = reinitialize_variables (newT,newV,newB)
%
% This function realloc global variables that describe the grid. 
%
% Input of this function are
%
% newT: the new size of variables that depend on nT
% newV: the new size of variables that depend on nV
% newB: the new size of variables that depend on nB
%

global TV TT B V VB TInfo BInfo BCircle
global nT nV nB

TV( nT+1:nT+newT , : ) = zeros ( newT , 3 );
TT( nT+1:nT+newT , : ) = zeros ( newT , 9 );
B( nB+1:nB+newB , : ) = zeros ( newB , 4 );
BInfo( nB+1:nB+newB , : ) = zeros ( newB , 3 );
V( nV+1:nV+newV , : ) = zeros ( newV , 2 );

% Struct for now are not initialized
%%% TInfo = ?;
%%% VB = ?;
%%% BCircle = ?:

return