% Sample 11

for i = 0: 99
    Domain.InputVertex(i+1,:) = [i/10 0];
end
Domain.InputVertex(101:106,:) = [ 10 0
    10 10
    6.2 10
    6 10
    5.8 10
    0 10];
Domain.Boundary.Values = 1:106 ;
Domain.Holes.Hole = [];
Domain.Segments.Segment = [];

BC.Values = 1;
BC.Boundary.Values = ones(1,106); 
BC.Holes.Hole = [];
BC.Segments.Segment = [];
BC.InputVertexValues = ones(1,106);

RefiningOptions.CheckArea = 'N';
RefiningOptions.CheckAngle = 'Y';
RefiningOptions.AreaValue = 0.1;   
RefiningOptions.AngleValue = 5;   
RefiningOptions.Subregions = [];

RefiningOptions.Subregions.Subregion(1).Type = 'rectangular';
RefiningOptions.Subregions.Subregion(1).Descriptors = [2 0 10 0];
RefiningOptions.Subregions.Subregion(1).CheckArea = 'N';
RefiningOptions.Subregions.Subregion(1).CheckAngle = 'Y';
RefiningOptions.Subregions.Subregion(1).AreaValue = []; 
RefiningOptions.Subregions.Subregion(1).AngleValue = 30; 

RefiningOptions.Subregions.Subregion(2).Type = 'elliptic';
RefiningOptions.Subregions.Subregion(2).Descriptors = [6 6 3 2];
RefiningOptions.Subregions.Subregion(2).CheckArea = 'Y';
RefiningOptions.Subregions.Subregion(2).CheckAngle = 'N';
RefiningOptions.Subregions.Subregion(2).AreaValue = 0.01; 
RefiningOptions.Subregions.Subregion(2).AngleValue = []; 

[geom] = btr30(Domain,BC,RefiningOptions);
BW_draw_grid (geom,1);