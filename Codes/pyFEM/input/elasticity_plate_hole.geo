//
// file elasticity_plate_hole.geo
//

// global dimensions

H=10;
V=10;
R=1.;

// mesh parameters
 
// fine
lc1=.5;
lc2=.05;

// coarse
//lc1=4;
//lc2=3.;

// point coordinates
Point(1) = {-H/2,-V/2,0,lc1};
Point(2) = {H/2,-V/2,0,lc1};
Point(3) = {H/2,V/2,0,lc1};
Point(4) = {-H/2,V/2,0,lc1};

Point(5) = {0,0,0,lc2};
Point(6) = {R,0,0,lc2};
Point(7) = {0,R,0,lc2};
Point(8) = {-R,0.,0,lc2};
Point(9) = {0,-R,0,lc2};

// lines and line loops
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(5) = {1,2,3,4};

Circle(6) = {6,5,7};
Circle(7) = {7,5,8};
Circle(8) = {8,5,9};
Circle(9) = {9,5,6};
Curve Loop(10) = {-6,-7,-8,-9};

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
Save 'elasticity_plate_hole.msh';
