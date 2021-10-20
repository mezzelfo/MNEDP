function [geom]=bbtr30(Domain,BC,RefiningOptions)
%
% _________________________________________________________________________
%
% BBTR - Barbera-Berrone TRiangler - version 3.0
%
% Politecnico di Torino
% 28 / 07 /08
% _________________________________________________________________________
% 
% This library can triangulate a plain surface that may contain every kind 
% of bond on the domain, as holes, concavities or inner segments.
%
% Parameters of the triangolation are the maximum area, the smallest angle  
% among the ones added by the library and the regions in which these
% parameters are active.
% 
% The output of this function is the set of data necessary to apply a 
% FiniteElement Method to the triangolation, and it includes boundary 
% conditions and element lists.
%
% Syntax and main features of BBTR 3.0 are exposed in an attached help pdf
% file
%


% -------------------------
% Load of global parameter
% -------------------------
global_parameter

% ----------------------------------------------------------
% First declaration and initialization of global variables
% ----------------------------------------------------------

global nT nV nB
global nKidnappedBorders KidnappedBorders EB nEB

global SplittingBordersLimits

% First declaration of number of elements
nT = 0;
nV = 0;
nB = 0;

% First declaration on Kidnapping List
nKidnappedBorders = 0;
KidnappedBorders = [];

% First declaration of Encroaching List
EB.Tree = [0 0];
EB.Br(1).Unsplitted = true;
EB.Br(1).Leaves = [];
EB.Br(1).nLeaves = 0;
EB.Br(1).KindOfSplitting = [];
EB.Br(1).SplittingThreshold = [];
EB.Br(1).SplittingLimits = [];
EB.nBr = 1;


% First declaration of every struct like a column array
VB = struct('n',[],'B',[],'V',[]);
VB(2,1).n = [];
TInfo = struct('Circumcenter',[],'Circumradius',[],'Area',[],'B',[]);
TInfo(2,1).Circumcenter = [];
BCircle = struct('Circumcenter',[],'r2',[]);
BCircle(2,1).Circumcenter = [];

% First extimation of sizes required to contain domain
[nTApprox nVApprox nBApprox] = first_size_variables_extimate ( length(Domain.InputVertex) );

% First declaration of every other variale
reinitialize_variables (nTApprox,nVApprox,nBApprox);




% --------------------------------------
% MAIN - Generation of triangulation
% --------------------------------------

% Modify SplittingBordersLimits to create fewer Subdomains;
TempSplittingBordersLimits = SplittingBordersLimits;
SplittingBordersLimits = Inf;

% Insert vertex in the Containing Box (CB)
[IntialVertexReference] = insert_vertex_in_box( Domain.InputVertex , size(Domain.InputVertex,1) );

% Add domain to the list of Borders to force into the domain
load_domain_in_kidnapping_list (Domain,BC,IntialVertexReference);

% Insert domain borders in triangulation
insert_domain_borders_in_box 

% Remove triangles out of domain (with exeption of particular ones)
[RemovedTriangles] = remove_box_triangles;

% Reset SplittingBordersLimits to chosen value;
SplittingBordersLimits = TempSplittingBordersLimits;

% Refining of the grid, using the refining options
refine_grid (RefiningOptions);

% Remove triangles in regions that had not to be triangulated (=outside triangles)
if RemovedTriangles == true
    
    % Modify SplittingBordersLimits to create fewer Subdomains;
    TempSplittingBordersLimits = SplittingBordersLimits;
    SplittingBordersLimits = Inf;

    % Remove triangles
    remove_outside_triangles
    
    % Reset SplittingBordersLimits to chosen value;
    SplittingBordersLimits = TempSplittingBordersLimits;
    
end

% Evaluate boundary conditions
[Di Ne Nodelist] = evaluate_boundary_conditions (Domain,BC);

% Transform arrays in simplier structures
[geom] = assign_variables_to_geom (Domain,BC,Di,Ne,Nodelist);

return