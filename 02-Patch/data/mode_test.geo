// 2D cantilever domain for Gmsh
// Shape:
// (0,0) -> (L,0) -> (L,H) -> (0,H)

SetFactory("Built-in");

// Geometry parameters
L  = 0.5;   // total length in x
lc = 0.8;  // target mesh size

// Points
Point(1) = {0, 0, 0, lc};
Point(2) = {3*L, 0, 0, lc};
Point(3) = {3*L, L, 0, lc};
Point(4) = {2*L, L, 0, lc};
Point(5) = {2*L, 2*L, 0, lc};
Point(6) = {L, 2*L, 0, lc};
Point(7) = {L, L, 0, lc};
Point(8) = {0, L, 0, lc};

// Boundary lines
Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 5}; // left
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

// Surface
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Physical groups for FEM/FEniCS
Physical Curve("bottom")         = {1}; // P1->P2: y=0,  base bottom
Physical Curve("right")          = {2}; // P2->P3: x=3L, outer right wall
Physical Curve("shoulder_right") = {3}; // P3->P4: y=L,  right ledge
Physical Curve("web_right")      = {4}; // P4->P5: x=2L, right side of web
Physical Curve("top")            = {5}; // P5->P6: y=2L, top of web
Physical Curve("web_left")       = {6}; // P6->P7: x=L,  left side of web
Physical Curve("shoulder_left")  = {7}; // P7->P8: y=L,  left ledge
Physical Curve("left")           = {8}; // P8->P1: x=0,  outer left wall
Physical Surface("domain") = {1};

// Mesh options
Mesh.Algorithm = 6; // Frontal-Delaunay for 2D
Mesh.ElementOrder = 2; // T6 triangles
Mesh.SecondOrderLinear = 1;
Mesh.MshFileVersion = 4.1;
