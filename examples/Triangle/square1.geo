h = 0.1;                         // Use a large value (we want the least possible number of elements)
Point(1) = {3, 0, 0, h};       
Point(2) = {3, 2, 0, h};
Point(3) = {0, 2, 0, h};
Point(4) = {0, 0, 0, h};
Line(12) = {1,2};  // box edges
Line(23) = {2,3};
Line(34) = {3,4};
Line(41) = {4,1};
//Line(42) = {4,2};  // diagonal line
//Curve Loop(1) = {12,23,34,41};  // the whole box
//Curve Loop(2) = {12,-42,41}; // the lower triangular region
Curve Loop(1234) = {12, 23, 34, 41};
//Plane Surface(19) = {1,2};  // upper triangular region = whole box - lower triangular region
//Plane Surface(20) = {2};  // lower triangular region
Plane Surface(111) = {1234};

Physical Surface("piastra",1) = {111}; 
Physical Point("Fixed") = {3,4};
Physical Point("Applied vertical force") = {2};
Physical Point("Horizontal roller") = {1};
Recombine Surface{:};

// To reproduce the mesh in two_triangles.py I define two separate plane 
// surfaces that are identical to the elements I'd like to obtain. 
// The two plane surfaces are defined by 1) the box minus the lower triangular 
// region and 2) the lower triangular region.
// GMSH doesn't seem to be able to build such a mesh without the 
// internal diagonal line (in that case, the mesh is symmetric with 4 T3 elements)

// NOTE: do no use tags to identify a physical group as this will not be understood by pyFEM
// - GOOD: Physical Surface("box") = {19,20}; 
// - BAD: Physical Surface("box",666) = {19,20}; and then use 666 in MaterialsSet in the input script


