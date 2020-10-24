h = 1;                         // Characteristic length of a mesh element
Point(1) = {0, 0, 0, 0.2*h};
Point(2) = {1, 0, 0, 0.2*h};
Point(3) = {1, 1, 0, 0.3*h};
Point(4) = {0, 1, 0, 0.4*h};
Point(5) = {0.350000, 0.100000, 0, 0.1*h};
Point(6) = {0.350000, 0.350000, 0, h};
Point(7) = {0.100000, 0.350000, 0, 0.1*h};
Point(8) = {0.85    , 0.85    , 0, h};
Line(1) = {1, 2}; // lower edge
Line(2) = {2, 3};
Line(3) = {3, 4}; // upper edge
Line(4) = {4, 1};
Circle(5) = {5,6,7};
Line(6) = {7,8};
Line(7) = {8,5};
Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5,6,7};

Plane Surface(8) = {1, 2}; // 8 and 9 are the elementary entity tags
Plane Surface(9) = {2};    // related to the two Line Loops

Physical Surface("matrix", 112) = {8}; 
Physical Surface("inclusion", 2) = {9}; 
Physical Curve("Dirichlet", 3) = {1}; 
Physical Curve("Neumann", 4) = {3};
Physical Point("Points",60) = {4,2,3};
//Physical Point("Dirichlet",60) = {4,2,3};

// Physical tags can be shown on the mesh by means of the GMSH GUI:
// 112 for the matrix (= domain - inclusion)
// 2 for the inclusion
// 3 for lower edge 
// 4 for upper edge
// 60 for three corner nodes
// 0 means that there is no tag.
// The physical tag values are arbitrary (and they have to be different from each other)
