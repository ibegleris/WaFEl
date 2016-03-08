cl__1 = 1e+22;
Point(1) = {-5e-05, -5e-05, -0, 1e+22};
Point(2) = {5e-05, 5e-05, -0, 1e+22};
Point(3) = {0, 0, -0, 1e+22};
Point(4) = {-5e-05, 5e-05, -0, 1e+22};
Point(5) = {5e-05, -5e-05, -0, 1e+22};
Point(7) = {1e-05, -0, -0, 1e+22};
Point(8) = {2e-05, -0, -0, 1e+22};
Line(1) = {4, 1};
Line(2) = {1, 5};
Line(3) = {5, 2};
Line(4) = {2, 4};
Circle(5) = {7, 3, 7};
Circle(6) = {8, 3, 8};

DefineConstant[ lc = { 0.1, Path "Gmsh/Parameters"}];
