// 2D Cook's membrane/domain for Gmsh
// Shape:
// (0,0) -> (L,t1) -> (L,t2) -> (0,t1)

SetFactory("Built-in");

// Geometry parameters
L  = 48;   // total length in x
t1 = 44;
t2 = 16;
t3 = t1 + t2;
lc = 0.5;  // target mesh size; smaller values create finer meshes

// Points
Point(1) = {0, 0, 0, lc};
Point(2) = {L, t1, 0, lc};
Point(3) = {L, t3, 0, lc};
Point(4) = {0, t1, 0, lc};

// Boundary lines
Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical groups 
Physical Curve("bottom")  = {1};
Physical Curve("right")   = {2};
Physical Curve("top") = {3};
Physical Curve("left") = {4};
Physical Surface("domain") = {1};

// Mesh options
Mesh.Algorithm = 6; // Frontal-Delaunay for 2D
Mesh.ElementOrder = 2; // Generate quadratic triangles (T6)
Mesh.MshFileVersion = 4.1;
