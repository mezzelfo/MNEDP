function []=global_parameter()
%
% []=global_parameter()
%
% This function sets the global parameters used in the library
%


global tolerance XBoxEnlarge YBoxEnlarge GenerationUntilResearch LimAngleToForceInGrid
global FirstExtimateTmultiplier FirstExtimateBmultiplier FirstExtimateVmultiplier
global FirstMaxRefiningRoof IncreasingRefiningRoof SplittingBordersLimits
global SafeMeshing RecursiveEncroachingLimit RefiningFolderName

% Tolerance used in the programme
tolerance = 1e-12;

% If this parameter is activated the program becomes slower but can
% triangulate also domains that otherwise would crash program (SafeMeshing
% must be activated if there are boundary borders very close to other ones)
SafeMeshing = true;

% Parameter used to create the Containing Box
XBoxEnlarge = 0.5;
YBoxEnlarge = 0.5;

% Generation until prioritary triangle research arrives to search
% (0 = Tstart, 1= neighbours, must be >1)
GenerationUntilResearch = 10;

% Triangle outside and adjacent to the domain with angle smaller than
% LimAngleToForceInGrid are forced in grid. It must be:
% 45 <= LimAngleToForceInGrid <= 90
LimAngleToForceInGrid = 45; %[°]

% Multipliers used to give a first extimate of dimension of elements of the
%grid
FirstExtimateTmultiplier = 8;
FirstExtimateBmultiplier = 8;
FirstExtimateVmultiplier = 6;

% Multiplier used to give a first extimate of Refing Roof
FirstMaxRefiningRoof = 1/5;

% Multiplier used to increase the Refining Roof after that it has been
% previously reached
IncreasingRefiningRoof = 3;

% Number of borders after which a domain subregion is splitted in two
SplittingBordersLimits = 10;

% Recursion stack limit
RecursionStackLimit = 1000;

% Recursion limit for small angles
RecursiveEncroachingLimit = 4;

% The fold that contain refining methods
RefiningFolderName = 'methods';



% -------------------------------------
% DO NOT MODIFY LINES BELOW THIS POINT
% -------------------------------------

set(0,'RecursionLimit',RecursionStackLimit)

return