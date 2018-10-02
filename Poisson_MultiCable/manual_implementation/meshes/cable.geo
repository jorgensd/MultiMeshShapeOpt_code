r = 100;
res = 4.0;
N = r*0.00025/res;  // Resolution
Point(1) = {0, 0, 0, 3*N};
Point(2) = {r*0.012, 0, 0, 3*N};
Point(3) = {r*-0.012, 0, 0, 3*N};
Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 3};
Physical Line(6) = {2, 1};
Line Loop(7) = {2, 1};
Plane Surface(8) = {7};
Physical Surface(10) = {8};
