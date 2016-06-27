a = DefineNumber[ 0.0002, Name "Parameters/a" ];
b = DefineNumber[ 0.0002, Name "Parameters/b" ];
rcore = DefineNumber[ 1e-5, Name "Parameters/rcore" ];
rclad = DefineNumber[ 10e-5, Name "Parameters/rclad" ];
lam = DefineNumber[ 1.55e-6, Name "Parameters/lam" ];
num = DefineNumber[  5, Name "Parameters/num" ];
sqrt22 = DefineNumber[  0.7071067811865476, Name "Parameters/sqrt22" ];

Point(1) = {0, 0, 0, lam/num};

Point(2) = {-a, -b, 0, 100*lam};
Point(3) = {-a, b, 0, 100*lam};
Point(4) = {a, b, 0, 100*lam};
Point(5) = {a, -b, 0, 100*lam};



Point(7) = {sqrt22*rcore, sqrt22*rcore, 0, lam/num};
Point(8) = {-sqrt22*rcore, sqrt22*rcore, 0,lam/num};
Point(9) = {-sqrt22*rcore, -sqrt22*rcore, 0, lam/num};
Point(10) ={sqrt22*rcore, -sqrt22*rcore, 0, lam/num};


Point(11) = {sqrt22*rclad, sqrt22*rclad, 0,50*lam};
Point(12) = {-sqrt22*rclad, sqrt22*rclad, 0, 50*lam};
Point(13) = {-sqrt22*rclad, -sqrt22*rclad, 0, 50*lam};
Point(14) = {sqrt22*rclad, -sqrt22*rclad, 0, 50*lam};


Circle(1) = {7, 1, 8};
Circle(2) = {8, 1, 9};
Circle(3) = {9, 1, 10};
Circle(4) = {10, 1, 7};


Circle(5) = {11, 1, 12};
Circle(6) = {12, 1, 13};
Circle(7) = {13, 1, 14};
Circle(8) = {14, 1, 11};

Line(9) = {2, 5};
Line(10) = {5, 4};
Line(11) = {4, 3};
Line(12) = {3, 2};


Line(13) = {2, 13};
Line(14) = {3, 12};
Line(15) = {4, 11};
Line(16) = {5, 14};
Line(17) = {13, 9};
Line(18) = {12, 8};
Line(19) = {11, 7};
Line(20) = {14, 10};
Line Loop(21) = {12, 13, -6, -14};
Plane Surface(22) = {21};
Line Loop(23) = {11, 14, -5, -15};
Plane Surface(24) = {23};
Line Loop(25) = {15, -8, -16, 10};
Plane Surface(26) = {25};
Line Loop(27) = {16, -7, -13, 9};
Plane Surface(28) = {27};
Line Loop(29) = {6, 17, -2, -18};
Plane Surface(30) = {29};
Line Loop(31) = {5, 18, -1, -19};
Plane Surface(32) = {31};
Line Loop(33) = {8, 19, -4, -20};
Plane Surface(34) = {33};
Line Loop(35) = {20, -3, -17, 7};
Plane Surface(36) = {35};
Line Loop(37) = {2, 3, 4, 1};
Plane Surface(38) = {37};

