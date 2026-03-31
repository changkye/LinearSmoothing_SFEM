// 2D cantilever domain for Gmsh
// Shape:
// (0,0) -> (L,0) -> (L,H) -> (0,H)

SetFactory("Built-in");

// Geometry parameters
L  = 1;   // total length in x
lc = 0.25;  // target mesh size

// Points
Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, L, 0, lc};
Point(4) = {0, L, 0, lc};

// Boundary lines
Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical groups for FEM/FEniCS
Physical Curve("bottom")         = {1}; // P1->P2: y=0,  base bottom
Physical Curve("right")          = {2}; // P2->P3: x=3L, outer right wall
Physical Curve("top")            = {3}; // P5->P6: y=2L, top of web
Physical Curve("left")           = {4}; // P8->P1: x=0,  outer left wall
Physical Surface("domain") = {1};

// Mesh options
Mesh.Algorithm = 6; // Frontal-Delaunay for 2D
Mesh.ElementOrder = 2; // T6 triangles
Mesh.SecondOrderLinear = 1;
Mesh.MshFileVersion = 4.1;
