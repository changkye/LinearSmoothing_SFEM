// 2D Cook's membrane for Gmsh
// Shape:
// (0,0) -> (48,44) -> (48,60) -> (0,44)

SetFactory("Built-in");

// Geometry parameters
L  = 48.0;  // right edge x-coordinate
H1 = 44.0;  // top-left y-coordinate
H2 = 60.0;  // top-right y-coordinate
lc = 1.0;   // target mesh size

// Points
Point(1) = {0, 0, 0, lc};
Point(2) = {L, 44, 0, lc};
Point(3) = {L, H2, 0, lc};
Point(4) = {0, H1, 0, lc};

// Boundary lines
Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical groups
Physical Curve("bottom") = {1};
Physical Curve("right") = {2};
Physical Curve("top") = {3};
Physical Curve("left") = {4};
Physical Surface("domain") = {1};

// Mesh options
Mesh.Algorithm = 6; // Frontal-Delaunay for 2D
Mesh.ElementOrder = 2; // T6 triangles
Mesh.SecondOrderLinear = 1;
Mesh.MshFileVersion = 4.1;
