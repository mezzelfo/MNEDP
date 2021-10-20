% Sample 12

Domain.InputVertex = [ 0 0
    5 0
    5 5
    0 5
    1 1
    4 1.25
    1 1.5
    4 4
    4 3.5
    1 3.75];
Domain.Boundary.Values = 1:4 ;
Domain.Holes.Hole(1).Values = 5:7;
Domain.Holes.Hole(2).Values = 8:10;
Domain.Segments.Segment(1).Values = [4 10];
Domain.Segments.Segment(2).Values = [10 6];
Domain.Segments.Segment(3).Values = [6 2];

BC.Values = [2 4 6 8 10 12];
BC.Boundary.Values = [1 2 1 2]; 
BC.Holes.Hole(1).Values = [3 3 3];
BC.Holes.Hole(2).Values = [4 4 4];
BC.Segments.Segment(1).Values = [5];
BC.Segments.Segment(2).Values = [6];
BC.Segments.Segment(3).Values = [7];
BC.InputVertexValues = [1 1 0 0 3 3 3 0 0 0];


RefiningOptions.CheckArea = 'Y';
RefiningOptions.CheckAngle = 'N';
RefiningOptions.AreaValue = 0.1;   
RefiningOptions.AngleValue = [];   
RefiningOptions.Subregions = [];


[geom] = btr30(Domain,BC,RefiningOptions);
BW_draw_grid (geom,1);