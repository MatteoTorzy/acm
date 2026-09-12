// global dimensions

H=10;
V=10;
R1=0.05;
R2=0.5;

// mesh parameters
lc1=.5;
lc2=.1;
lc3=.001;

// point coordinates
Point(1) = {-H/2,-V/2,0,lc1};
Point(2) = {H/2,-V/2,0,lc1};
Point(3) = {H/2,V/2,0,lc1};
Point(4) = {-H/2,V/2,0,lc1};

Point(5) = {0,0,0,lc2};
Point(6) = {R1,0,0,lc2};
Point(7) = {0,R2,0,lc3};
Point(8) = {-R1,0.,0,lc2};
Point(9) = {0,-R2,0,lc3};

// lines and line loops
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

Ellipse(6) = {6,5,7,7};
Ellipse(7) = {7,5,7,8};
Ellipse(8) = {8,5,7,9};
Ellipse(9) = {9,5,7,6};
Line Loop(10) = {-6,-7,-8,-9};

// surface
Plane Surface(1) = {5,10};

// physical entities
Physical Line("RightEdge") = {2}; 
Physical Line("LeftEdge") = {4}; 
Physical Surface("Plate") = {1}; 

// creates second order mesh and saves
Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'elasticity_plate_ellipse.msh';
