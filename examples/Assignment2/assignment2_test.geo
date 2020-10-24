// Mesh size
h = 3.5;

// Constants
rd = 2.5;     // Hole radius.
lx = 20;    // Plate width.
ly = 20;    // Plate height.

// Points definition
Point(1) = {     0,      0, 0, h};
Point(2) = {    rd,      0, 0, h};
Point(3) = {lx*0.5,      0, 0, h};
Point(4) = {lx*0.5, ly*0.5, 0, h};
Point(5) = {     0, ly*0.5, 0, h};
Point(6) = {     0,     rd, 0, h};

// Line definition
Circle(62) = {6,1,2};   //62 to have counter clock-wise direction
Line(23)   = {2,3};
Line(34)   = {3,4};
Line(45)   = {4,5};
Line(56)   = {5,6};

// Line loop
Line Loop(1) = {62, 23, 34, 45, 56};

// Surface generation
Plane Surface(1) = {1};

// Physical surface
Physical Surface("piastra", 100) = {1};

// Physical lines
//Physical Curve("x-symmetry",        101) = {23}; 
//Physical Curve("y-symmetry",        102) = {56}; 
Physical Curve("EnforcedDisplacement",   103) = {23, 34, 45, 56, 62};
