L = 0.6;

// Points
Point(1) = {0, 0, 0, L};
Point(2) = {1, 0, 0, L};
Point(3) = {0, 1, 0, L};

// Lines
Line(1) = {3, 1};
Line(2) = {1, 2};
Circle(3) = {2, 1, 3};

// Surfaces
Line Loop(1) = {2, 3, 1};
Plane Surface(1) = {1};

// Physical groups
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Point(3) = {3};
Physical Surface(4) = {1};
