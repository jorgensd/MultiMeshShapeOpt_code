r = 100; //Scaling
res = 4.0;
r0 = r*0.002;
r1 = r*0.00255; // Rubber radius
r2 = r*0.003;   // Halo radius
N = r*0.00025/res;  // Resolution

Point(1) = {0, 0, 0, N};
Point(2) = { r0, 0, 0, N};
Point(3) = {-r0, 0, 0, N};
Point(4) = { r1, 0, 0, N};
Point(5) = {-r1, 0, 0, N};
Point(6) = { r2, 0, 0, 3*N};
Point(7) = {-r2, 0, 0, 3*N};

Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 3};
Circle(3) = {5, 1, 4};
Circle(4) = {4, 1, 5};
Circle(5) = {7, 1, 6};
Circle(6) = {6, 1, 7};
Line Loop(7) = {6, 5};
Line Loop(8) = {4, 3};
Plane Surface(9) = {7, 8};
Line Loop(10) = {2, 1};
Plane Surface(11) = {8, 10};
Plane Surface(12) = {10};
Physical Surface(13) = {9};
Physical Surface(14) = {11};
Physical Surface(15) = {12};
Physical Line(16) = {2, 1};
Physical Line(17) = {3, 4};
Physical Line(18) = {5, 6};