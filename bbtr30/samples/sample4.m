% Sample 11

Domain.InputVertex = [ 0 0
    5 0
    5 5
    0 5];

Domain.Boundary.Values = [1 2 3 4] ;
Domain.Holes.Hole = [];
Domain.Segments.Segment = [];

BC.Values = 1;
BC.Boundary.Values = [1 1 1 1]; 
BC.Holes.Hole = [];
BC.Segments.Segment = [];
BC.InputVertexValues = [1 1 1 1];

RefiningOptions.CheckArea = 'Y';
RefiningOptions.CheckAngle = 'N';
RefiningOptions.AreaValue = 0.1;   
RefiningOptions.AngleValue = 30;   
RefiningOptions.Subregions = [];

% RefiningOptions.Subregions.Subregion(1).Type = 'rectangular';
% RefiningOptions.Subregions.Subregion(1).Descriptors = [4 1 3 2];
% RefiningOptions.Subregions.Subregion(1).CheckArea = 'Y';
% RefiningOptions.Subregions.Subregion(1).CheckAngle = 'N';
% RefiningOptions.Subregions.Subregion(1).AreaValue = 0.01; 
% RefiningOptions.Subregions.Subregion(1).AngleValue = []; 


[geom] = btr30(Domain,BC,RefiningOptions);
draw_grid (geom,1);