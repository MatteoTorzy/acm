lc = .5;
rad1=2.;
rad2=3.;
ver=5;

// point coordinates
Point(1) = {0.,ver,0.,lc};
Point(2) = {rad2,ver,0.,lc};
Point(3) = {0.,ver+rad2,0.,lc};
Point(4) = {-rad2,ver,0.,lc};
Point(5) = {0.,ver-rad2,0,lc};
Point(6) = {rad1,ver,0.,lc};
Point(7) = {0.,ver+rad1,0.,lc};
Point(8) = {-rad1,ver,0.,lc};
Point(9) = {0.,ver-rad1,0,lc};
Point(10) = {-rad2,rad2/2,lc};
Point(11) = {rad2,-rad2/2,0,lc};

// lines and circles
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};
Line(9) = {10,11};

Line Loop(10) = {1,2,3,4};
Line Loop(11) = {5,6,7,8};

// surfaces
Plane Surface(1) = {10,11};

Physical Surface("Tyre") = {1}; 
Physical Line("Border") = {1,2,3,4}; 
Physical Line("Contact") = {9}; 

// creates second order mesh and saves

Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'dynamics_tyreincl.msh';
