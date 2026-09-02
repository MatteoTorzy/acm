//SetFactory("OpenCASCADE");

// -------------------------
// Parameters
// -------------------------
Lx = 6.0;
Ly = 4.0;
xm = Lx/3;

xc = 4.;   // hole center x (right half)
yc = 2.0;   // hole center y
r  = 0.6;   // hole radius

lc1 = 0.4;
lc2 = 0.05;
lc3 = 0.1;

// -------------------------
// Geometry points
// -------------------------
Point(1) = {0,  0,  0, lc1};
Point(2) = {xm, 0,  0, lc3};
Point(3) = {Lx, 0,  0, lc1};
Point(4) = {Lx, Ly, 0, lc1};
Point(5) = {xm, Ly, 0, lc3};
Point(6) = {0,  Ly, 0, lc1};

// Hole points
Point(10) = {xc,    yc,    0, lc2};
Point(11) = {xc+r,  yc,    0, lc2};
Point(12) = {xc,    yc+r,  0, lc2};
Point(13) = {xc-r,  yc,    0, lc2};
Point(14) = {xc,    yc-r,  0, lc2};

// -------------------------
// Outer boundary lines
// -------------------------
Line(1) = {1, 2};   // bottom-left
Line(2) = {2, 3};   // bottom-right
Line(3) = {3, 4};   // right
Line(4) = {4, 5};   // top-right
Line(5) = {5, 6};   // top-left
Line(6) = {6, 1};   // left
Line(7) = {2, 5};   // middle vertical split

// -------------------------
// Hole circle
// -------------------------
Circle(8)  = {11, 10, 12};
Circle(9)  = {12, 10, 13};
Circle(10) = {13, 10, 14};
Circle(11) = {14, 10, 11};

// -------------------------
// Left subdomain
// -------------------------
Curve Loop(20) = {1, 7, 5, 6};
Plane Surface(30) = {20};

// -------------------------
// Right subdomain with hole
// -------------------------
Curve Loop(21) = {2, 3, 4, -7};
Curve Loop(22) = {-8, -9, -10, -11};
Plane Surface(31) = {21, 22};

// -------------------------
// Physical groups
// -------------------------
Physical Surface("LeftPart")  = {30};
Physical Surface("RightPart") = {31};
Physical Curve("LeftEdge")   = {6};
Physical Curve("RightEdge")  = {3};
Physical Curve("LowerEdge")   = {1,2};
Physical Curve("UpperEdge")  = {4,5};
Physical Curve("Hole")       = {-8, -9, -10, -11};

// creates second order mesh and saves
Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'example.msh';
