% Sample 12

Domain.InputVertex = [ 0 0
    1 0
    1 1
    0 1];
Domain.Boundary.Values = [1 2 3 4] ;
Domain.Holes.Hole = [];
Domain.Segments.Segment = [];


BC.Boundary.Values = [1 2 4 3]; 
BC.InputVertexValues = [5 7 0 1];
BC.Values = [0.0 20.0 0.12 0.0 0.654 -1 0.33]
BC.Holes.Hole = [];
BC.Segments.Segment = [];

RefiningOptions.CheckArea = 'Y';
RefiningOptions.CheckAngle = 'N';
RefiningOptions.AreaValue = 0.2;   
RefiningOptions.AngleValue = [];   
RefiningOptions.Subregions = [];


[geom] = btr30(Domain,BC,RefiningOptions);
draw_grid (geom,1);