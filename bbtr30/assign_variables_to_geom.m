function [geom] = assign_variables_to_geom(Domain,BC,Di,Ne,Nodelist);
%
% [geom] = assign_variables_to_geom(Domain,BC,Di,Ne,Nodelist);
%
% This function assign to geom struct the variables involved in the grid
%
% Input & Output of this function are:
%
% Domain:
%   Domain.Boundary.Values: references to boundary verteces
%   Domain.Hole.Values: references to verteces belonging to a hole
%   Domain.Segment.Values: references to verteces belonging to a segment
%   Domain.InputVertex coordinates of input verteces
%   (*see attachments to find usage informations on Domain variable)
%
% BC:
%   BC.Boundary.Values: references to boundary BC
%   BC.Hole.Values: references to holes BC
%   BC.Segment.Values: references to segments BC
%   BC.Values: references to BC values
%   BC.InputVertexValues: references to BC to be imposed in InputVertex
%   (*see attachments to find usage informations on BC variable)
%
% Nodelist: an array that keep information of BCs to impose in every vertex
%   Nodelist(i) = 0 : the i-vertex has no BC related to it
%   Nodelist(i) = j : the i-vertex has BC = j (where j is a reference to 
%       BC.Values). This condition must be a Dirichlet one
%
% Ne: this is a list of borders whose assigned BC are Neumann ones.
%   The i-esim row of Ne is [i,j], where i is the i-esim border and j is
%   the reference to the condition related to it.
%
% Di: this is a list of verteces whose assigned BC are Dirichlet ones.
%   The i-esim row of De is [i,j], where i is the i-esim vertex and j is
%   the reference to the condition related to it.
%
% geom: a structure containing every information on the triangulation
%       (*see attachments for further informations)

global TV TT B V VB TInfo BInfo BCircle
global nT nV nB
global EB

geom.elements.triangles = TV;
clear TV
geom.elements.coordinates = V;
clear V
geom.elements.borders = B;
clear B
geom.elements.neighbourhood = TT;
clear TT
geom.elements.vertexesneighbourhood = VB;
clear VB
geom.nelements.nTriangles = nT;
geom.nelements.nBorders = nB;
geom.nelements.nVertexes = nV;
geom.pivot.nodelist = Nodelist;
geom.pivot.Di = Di;
geom.pivot.Ne = Ne;
geom.support.TInfo = TInfo;
clear TInfo
geom.support.BInfo = BInfo;
clear BInfo
geom.support.BCircle = BCircle;
clear BCircle
geom.support.EB = EB;
clear EB
geom.input.Domain = Domain;
geom.input.BC = BC;


return
