// Gmsh project created on Thu Feb  2 06:13:15 2017
Point(1) = {-10, 10, 0, 1};
Point(2) = {2, 10, 0, 1};
Point(3) = {2, 0, 0, 1};
Point(4) = {0, 0, 0, 1};
Point(5) = {-10, 0, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
//+
Line Loop(6) = {1, 2, 3, 4, 5};
//+
Plane Surface(7) = {6};
//+
Physical Line("Freestream") = {5, 1};
//+
Physical Line("Outflow") = {2};
//+
Physical Line("Symmetric") = {4};
//+
Physical Line("Wall") = {3};
//+
Physical Surface("Flow") = {7};
//+
Field[1] = BoundaryLayer;
//+
Field[1].NodesList = {4, 3};
//+
Field[1].EdgesList = {3};
//+
Field[1].Quads = 1;
//+
Field[1].AnisoMax = 0;
//+
Field[1].thickness = 1.5;
//+
Field[1].FanNodesList = {4};
//+
Field[1].hwall_n = 0.01;
//+
Field[1].hfar = 1;
//+
Field[1].ratio = 1.05;
//+
Background Field = 1;
