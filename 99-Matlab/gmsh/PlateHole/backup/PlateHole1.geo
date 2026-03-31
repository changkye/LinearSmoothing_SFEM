h  = 1.0;
lx = 1;
ly = 1;
r  = 0.5;

Point(1)  = {0,0,0,h};
Point(2)  = {r,0,0,h};
Point(3)  = {0,r,0,h};
Point(4)  = {-r,0,0,h};
Point(5)  = {0,-r,0,h};
Point(6)  = {lx,0,0,h};
Point(7)  = {lx,ly,0,h};
Point(8)  = {0,ly,0,h};
Point(9)  = {-lx,ly,0,h};
Point(10) = {-lx,0,0,h};
Point(11) = {-lx,-ly,0,h};
Point(12) = {0,-ly,0,h};
Point(13) = {lx,-ly,0,h};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Line(5)   = {6,7};
Line(6)   = {7,8};
Line(7)   = {8,9};
Line(8)   = {9,10};
Line(9)   = {10,11};
Line(10)  = {11,12};
Line(11)  = {12,13};
Line(12)  = {13,6}; 

Line Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12};

Plane Surface(1) = {1};
Physical Surface(1) = {1};

Physical Line(13) = {10, 11, 12, 5, 8, 9}; \\--- BC
Physical Line(14) = {7, 6}; \\--- load
