h  = 0.5;
lx = 2;
ly = 2;
r  = 1;

Point(1) = {0,0,0,h};
Point(2) = {r,0,0,h};
Point(3) = {lx,0,0,h};
Point(4) = {lx,ly,0,h};
Point(5) = {0,ly,0,h};
Point(6) = {0,r,0,h};

Circle(1) = {6,1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};

Line Loop(1) = {1,2,3,4,5};

Plane Surface(1) = {1};
Physical Surface(1) = {1};
