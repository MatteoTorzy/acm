H=1;
V=20;

lc1=.3;
lc2=.3;

Point(1) = {0,0,0,lc1};
Point(2) = {H,0,0,lc1};
Point(3) = {H,V/2-.1,0,lc2};
Point(4) = {H,V,0,lc1};
Point(5) = {0,V,0,lc1};
Point(6) = {0,V/2,0,lc2};
Point(7) = {H,V/2+.1,0,lc2};

Line (1) = {1,2};
Line (2) = {2,3};
Line (3) = {3,7};
Line (4) = {7,4};
Line (5) = {4,5};
Line (6) = {5,6};
Line (7) = {6,1};

Line Loop (8) = {1,2,3,4,5,6,7};

Plane Surface (1) = {8};

Physical Line ("Lower") = {1};
Physical Line ("Upper") = {5};
Physical Line ("Load") = {3};
Physical Surface ("Solid") = {1};

// creates second order mesh and saves
Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'buckling.msh';

