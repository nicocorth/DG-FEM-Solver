Point(1) = {-1,-1,-1};
Point(2) = {-1, 1,-1};
Point(3) = { 1, 1,-1};
Point(4) = { 1,-1,-1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Recombine Surface{1};
Extrude {0, 0, 2} {
  Surface{1}; Layers{5}; Recombine;
}
//+
Transfinite Curve {2, 1, 4, 3} = 5 Using Progression 1;
//+
Physical Surface("Sinusoidal") = {25};
//+
Physical Surface("Wall") = {26, 13, 1, 21};
//+
Physical Surface("Opening") = {17};
//+
Physical Volume("Domain") = {1};
