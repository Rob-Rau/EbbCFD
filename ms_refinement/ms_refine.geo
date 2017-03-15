// Gmsh project created on Thu Feb  2 06:13:15 2017
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//+
Line Loop(6) = {1, 2, 3, 4};
//+
Plane Surface(7) = {6};
//+
Physical Line("Freestream") = {1, 2, 3, 4};
//+
Physical Surface("Flow") = {7};

