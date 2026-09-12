
// global dimensions
H=10;
V=0.5;

// mesh parameters

lc1=V/4;

// point coordinates
Point(1) = {0,0,0,lc1};
Point(2) = {H,0,0,lc1};
Point(3) = {H,V,0,lc1};
Point(4) = {0,V,0,lc1};

// lines and line loops
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

// surface
Plane Surface(1) = {5};

// physical entities
Physical Line("LeftEdge") = {4};   // left 
Physical Line("RightEdge") = {2};   // right
Physical Surface("Tyre") = {1};  // surface

// creates second order mesh and saves
Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'modal_beamCC.msh';
