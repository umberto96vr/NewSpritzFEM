// Gmsh project created on Fri May 01 09:59:27 2020
h = 1;
L = 10;

Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {0.5*L, 0.8660254038*L, 0, h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};

Curve Loop(1) = {1,2,3};
Plane Surface(1) = {1};

Physical Surface(1) = {1};