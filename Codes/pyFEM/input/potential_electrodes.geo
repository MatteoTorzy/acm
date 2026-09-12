// Parametric capacitor in air

lc = 0.2;

// Electrode parameters
electrode_length    = 4.0;
electrode_thickness = 0.4;
electrode_gap       = 2.0;

// Air margins
air_margin_parallel = 3.0;
air_margin_gap      = 3.0;

// Derived dimensions
Lx = electrode_gap + 2*electrode_thickness + 2*air_margin_gap;
Ly = electrode_length + 2*air_margin_parallel;

// Electrode positions
xL = air_margin_gap;
xR = air_margin_gap + electrode_thickness + electrode_gap;

y0 = air_margin_parallel;
y1 = air_margin_parallel + electrode_length;

// Outer air rectangle
Point(1) = {0,0,0,lc};
Point(2) = {Lx,0,0,lc};
Point(3) = {Lx,Ly,0,lc};
Point(4) = {0,Ly,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};

// Left electrode hole
Point(11) = {xL,y0,0,lc};
Point(12) = {xL+electrode_thickness,y0,0,lc};
Point(13) = {xL+electrode_thickness,y1,0,lc};
Point(14) = {xL,y1,0,lc};

Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,11};
Curve Loop(2) = {11,12,13,14};

// Right electrode hole
Point(21) = {xR,y0,0,lc};
Point(22) = {xR+electrode_thickness,y0,0,lc};
Point(23) = {xR+electrode_thickness,y1,0,lc};
Point(24) = {xR,y1,0,lc};

Line(21) = {21,22};
Line(22) = {22,23};
Line(23) = {23,24};
Line(24) = {24,21};
Curve Loop(3) = {21,22,23,24};

// Air domain
Plane Surface(1) = {1,2,3};

// Physical groups
Physical Surface("Air") = {1};

Physical Curve("OuterBoundary") = {1,2,3,4};
Physical Curve("Electrode1") = {11,12,13,14};
Physical Curve("Electrode2") = {21,22,23,24};


// Second-order triangular mesh in Gmsh msh v2 format.
Mesh.MshFileVersion = 2;
Mesh.ElementOrder = 2;
Mesh 2;
Save "potential_electrodes.msh";

