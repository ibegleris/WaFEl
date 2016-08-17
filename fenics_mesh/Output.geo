a = DefineNumber[ 0.0001, Name "Parameters/a" ];
b = DefineNumber[ 0.0001, Name "Parameters/b" ];
rcore = DefineNumber[ 2e-06, Name "Parameters/rcore" ];
rclad = DefineNumber[ 5e-05, Name "Parameters/rclad" ];
numlim = DefineNumber[ 1, Name "Parameters/numlim" ];
lam = DefineNumber[ 1.55e-06, Name "Parameters/lam" ];
elmin = DefineNumber[ numlim*lam, Name "Parameters/elmin" ];
sqrt2 = DefineNumber[ 0.7071067811865476, Name "Parameters/sqrt2" ];
elminclad = DefineNumber[ 70*numlim*lam, Name "Parameters/elminclad" ];
elminboundary = DefineNumber[ 100*numlim*lam, Name "Parameters/elminboundary" ];



Point(1) = {0, 0, 0, elmin};
Point(2) = {rcore, 0, 0,elmin};
Point(3) = {sqrt2*rcore, sqrt2*rcore, 0,elmin};
Point(4) = {0, rcore, 0,elmin};
Point(5) = {-sqrt2*rcore, sqrt2*rcore, 0,elmin};
Point(6) = {-rcore, 0, 0,elmin};
Point(7) = {-sqrt2*rcore, -sqrt2*rcore, 0,elmin};
Point(8) = {0, -rcore, 0,elmin};
Point(9) = {sqrt2*rcore, -sqrt2*rcore, 0,elmin};


Point(10) = {rclad, 0, 0,elminclad};
Point(11) = {sqrt2*rclad, sqrt2*rclad, 0,elminclad};
Point(12) = {0, rclad, 0,elminclad};
Point(13) = {-sqrt2*rclad, sqrt2*rclad, 0,elminclad};
Point(14) = {-rclad, 0, 0,elminclad};
Point(15) = {-sqrt2*rclad, -sqrt2*rclad, 0,elminclad};
Point(16) = {0, -rclad, 0,elminclad};
Point(17) = {sqrt2*rclad, -sqrt2*rclad, 0,elminclad};



Point(18) = {0, -b, 0,elminboundary};
Point(19) = {a, 0, 0,elminboundary};
Point(20) = {0, b, 0,elminboundary};
Point(21) = {-a, 0, 0,elminboundary};


Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 6};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 2};



Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 14};
Circle(13) = {14, 1, 15};
Circle(14) = {15, 1, 16};
Circle(15) = {16, 1, 17};
Circle(16) = {17, 1, 10};



//Line(17) = {18, 19};
//Line(18) = {19, 20};
//Line(19) = {20, 21};
//Line(20) = {21, 18};
Circle(17) = {19, 1, 20};
Circle(18) = {20, 1, 21};
Circle(19) = {21, 1, 18};
Circle(20) = {18, 1, 19};


Line Loop(21) = {20, 17, 18, 19};
Line Loop(22) = {12, 13, 14, 15, 16, 9, 10, 11};
Plane Surface(23) = {21, 22};
Line Loop(24) = {3, 4, 5, 6, 7, 8, 1, 2};
Plane Surface(25) = {22, 24};
Plane Surface(26) = {24};
//Transfinite Surface {23};
//Transfinite	Surface {25};
//Trnasfinite Surface {26};
//Recombine Surface {23};
//Recombine Surface {25};
//Recombine Surface {26};

