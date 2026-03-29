// 2D L-shaped beam/domain for Gmsh
// Shape:
// (0,0) -> (L,0) -> (L,t) -> (t,t) -> (t,H) -> (0,H)

SetFactory("Built-in");

// Geometry parameters
L  = 1.0;   // total length in x
H  = 1.0;   // total height in y
t  = 0.4;   // thickness of the L-arm
lc = 0.01;  // target mesh size

// Points
Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, t, 0, lc};
Point(4) = {t, t, 0, lc};
Point(5) = {t, H, 0, lc};
Point(6) = {0, H, 0, lc};

// Boundary lines
Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // inner horizontal
Line(4) = {4, 5}; // inner vertical
Line(5) = {5, 6}; // top
Line(6) = {6, 1}; // left

// Surface
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// Physical groups for FEM/FEniCS
Physical Curve("bottom")  = {1};
Physical Curve("right")   = {2};
Physical Curve("inner_h") = {3};
Physical Curve("inner_v") = {4};
Physical Curve("top")     = {5};
Physical Curve("left")    = {6};
Physical Surface("domain") = {1};

// Mesh options
Mesh.Algorithm = 6; // Frontal-Delaunay for 2D
Mesh.ElementOrder = 2; // T6 triangles
Mesh.SecondOrderLinear = 1;
Mesh.MshFileVersion = 4.1;
