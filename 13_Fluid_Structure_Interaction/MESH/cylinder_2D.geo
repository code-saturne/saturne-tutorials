// Set refinement
raf = 4.0;

// Add cylinder points
Radius = 0.8;
CenterX = 5.0;
CenterZ = 0.0;
CenterY = 0.0;
Point(1) = {CenterX, CenterY, CenterZ, 1.0};
Point(2) = {CenterX+Radius*Cos(Pi/4.), CenterY, CenterZ+Radius*Sin(Pi/4.), 1.0};
Point(3) = {CenterX-Radius*Cos(Pi/4.), CenterY, CenterZ+Radius*Sin(Pi/4.), 1.0};
Point(4) = {CenterX-Radius*Cos(Pi/4.), CenterY, CenterZ-Radius*Sin(Pi/4.), 1.0};
Point(5) = {CenterX+Radius*Cos(Pi/4.), CenterY, CenterZ-Radius*Sin(Pi/4.), 1.0};

// Add extrusion around cylinder
Radius_ext = 4.0;
Point(6) = {CenterX+Radius_ext*Cos(Pi/4.), CenterY, CenterZ+Radius_ext*Sin(Pi/4.), 1.0};
Point(7) = {CenterX-Radius_ext*Cos(Pi/4.), CenterY, CenterZ+Radius_ext*Sin(Pi/4.), 1.0};
Point(8) = {CenterX-Radius_ext*Cos(Pi/4.), CenterY, CenterZ-Radius_ext*Sin(Pi/4.), 1.0};
Point(9) = {CenterX+Radius_ext*Cos(Pi/4.), CenterY, CenterZ-Radius_ext*Sin(Pi/4.), 1.0};

// Add boundaries points
LenghtX = 38.0;
LenghtZ = 12.0;
LenghtY = 0.5;
Point(10) = {0., CenterY, CenterZ+Radius_ext*Sin(Pi/4.), 1.0};
Point(11) = {0., CenterY, CenterZ-Radius_ext*Sin(Pi/4.), 1.0};
Point(12) = {CenterX-Radius_ext*Cos(Pi/4.), CenterY, LenghtZ/2., 1.0};
Point(13) = {CenterX+Radius_ext*Cos(Pi/4.), CenterY, LenghtZ/2., 1.0};
Point(14) = {CenterX-Radius_ext*Cos(Pi/4.), CenterY, -LenghtZ/2., 1.0};
Point(15) = {CenterX+Radius_ext*Cos(Pi/4.), CenterY, -LenghtZ/2., 1.0};
Point(16) = {LenghtX, CenterY, CenterZ+Radius_ext*Sin(Pi/4.), 1.0};
Point(17) = {LenghtX, CenterY, CenterZ-Radius_ext*Sin(Pi/4.), 1.0};
Point(18) = {LenghtX, CenterY, LenghtZ/2., 1.0};
Point(19) = {LenghtX, CenterY, -LenghtZ/2., 1.0};
Point(20) = {0, CenterY, LenghtZ/2., 1.0};
Point(21) = {0, CenterY, -LenghtZ/2., 1.0};

// Add cylinder lines
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Add extrusion lines around cylinder
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
Line(9) = {3, 7};
Line(10) = {4, 8};
Line(11) = {5, 9};
Line(12) = {2, 6};

// Add boundary lines
Line(13) = {7, 12};
Line(14) = {6, 13};
Line(15) = {6, 16};
Line(16) = {9, 17};
Line(17) = {9, 15};
Line(18) = {8, 14};
Line(19) = {8, 11};
Line(20) = {7, 10};
Line(21) = {12, 20};
Line(22) = {14, 21};
Line(23) = {11, 21};
Line(24) = {10, 20};
Line(25) = {16, 18};
Line(26) = {17, 19};
Line(27) = {13, 18};
Line(28) = {15, 19};
Line(29) = {14, 15};
Line(30) = {17, 16};
Line(31) = {13, 12};
Line(32) = {10, 11};

// Set lines refinement around cylinder
Resol_cyl = 10.*raf;
Conv = Resol_cyl/(Pi*Radius);
Resol_ext = (Radius_ext-Radius)*Conv;
Transfinite Line {1, 2, 3, 4, 5, 6, 7, 8} = Resol_cyl Using Progression 1;
Transfinite Line {9, 10, 11, 12} = Resol_ext Using Progression 1.;

// Set lines refinement around boundaries
Resol_1 = (CenterX-Radius_ext*Cos(Pi/4.))*Conv*0.5;
Resol_2 = (LenghtZ/2.-(CenterZ+Radius_ext*Sin(Pi/4.)))*Conv*0.5;
Resol_3 = (LenghtX-(CenterX+Radius_ext*Cos(Pi/4.)))*Conv*0.5;
Transfinite Line {19, 20, 21, 22} = Resol_1 Using Progression 1;
Transfinite Line {13, 14, 24, 25, 23, 18, 17, 26} = Resol_2 Using Progression 1;
Transfinite Line {29, 30, 31, 32} = Resol_cyl Using Progression 1;
Transfinite Line {15, 16, 27, 28} = Resol_3 Using Progression 1;

// Define surfaces around cylinder
Line Loop(1) = {5, -9, -1, 12}; Plane Surface(1) = {1};
Line Loop(2) = {9, 6, -10, -2}; Plane Surface(2) = {2};
Line Loop(3) = {3, 11, -7, -10}; Plane Surface(3) = {3};
Line Loop(4) = {12, -8, -11, 4}; Plane Surface(4) = {4};
Transfinite Surface {1, 2, 3, 4};
Recombine   Surface {1, 2, 3, 4};

// Define boundary surfaces
Line Loop(5) = {21, -24, -20, 13}; Plane Surface(5) = {5};
Line Loop(6) = {23, -22, -18, 19}; Plane Surface(6) = {6};
Line Loop(7) = {32, -19, -6, 20}; Plane Surface(7) = {7};
Line Loop(8) = {7, 17, -29, -18}; Plane Surface(8) = {8};
Line Loop(9) = {31, -13, -5, 14}; Plane Surface(9) = {9};
Line Loop(10) = {27, -25, -15, 14}; Plane Surface(10) = {10};
Line Loop(11) = {16, 26, -28, -17}; Plane Surface(11) = {11};
Line Loop(12) = {15, -30, -16, 8}; Plane Surface(12) = {12};
Transfinite Surface {5, 6, 7, 8, 9, 10, 11, 12};
Recombine   Surface {5, 6, 7, 8, 9, 10, 11, 12};

// 2D-Sym Extrusion
Extr = 0.5;
Extrude {0, Extr, 0} {
  Surface{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Layers{1};
  Recombine;
}

// Inlet
Physical Surface(1) = {133, 173, 151};
// Outlet
Physical Surface(2) = {243, 287, 265};
// Cylinder walls
Physical Surface(3) = {119, 85, 49, 75};
// 2D-sym
Physical Surface(4) = {6, 8, 11, 12, 10, 9, 5, 7, 2, 3, 4, 1, 252, 296, 274, 208, 164, 186, 142, 230, 54, 120, 98, 76};
// Up and down
Physical Surface(5) = {129, 217, 239, 269, 203, 155};

Physical Volume(10) = {5, 7, 9, 2, 1, 6, 3, 4, 8, 10, 12, 11};
