Point(1) = {-1,-1,-1};
Point(2) = {-1, 1,-1};
Point(3) = { 1, 1,-1};
Point(4) = { 1,-1,-1};

Point(5) = {-1,-1, 1};
Point(6) = {-1, 1, 1};
Point(7) = { 1, 1, 1};
Point(8) = { 1,-1, 1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};//+
Curve Loop(1) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, -5, -9, 1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, -9, -4, 12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, -12, -3, 11};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {11, -6, -10, 2};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {1, 2, 3, 4};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {3, 1, 2, 5, 4, 6};
//+
Volume(1) = {1};
//+
Transfinite Curve {2, 10, 6, 11, 5, 1, 3, 7, 12, 4, 9, 8} = 3 Using Progression 1;
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
//+
Transfinite Surface {2};
//+
Transfinite Surface {4};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
Physical Volume("Domain") = {1};
//+
Physical Surface("Sinusoidal 1") = {2};
//+
Physical Surface("Wall 1") = {1};
//+
Physical Surface("Wall 2") = {3};
//+
Physical Surface("Wall 3") = {6};
//+
Physical Surface("Wall 4") = {5};
//+
Physical Surface("Opening") = {4};
