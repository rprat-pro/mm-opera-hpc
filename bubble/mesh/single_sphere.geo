size_rev=10e-6; // the size of the REV
mesh_size= size_rev*0.5e-1; 
scaling_mesh_bubble = 0.5e-1;
radius_bubble = 400e-9;

//+
Point(1) = {0, 0, 0, mesh_size*scaling_mesh_bubble};
Point(2) = {radius_bubble, 0, 0, mesh_size*scaling_mesh_bubble};
Point(3) = {-radius_bubble, 0, 0, mesh_size*scaling_mesh_bubble};
Point(4) = {0, 0, radius_bubble, mesh_size*scaling_mesh_bubble};
Point(5) = {0, 0, -radius_bubble, mesh_size*scaling_mesh_bubble};
Point(6) = {0, radius_bubble, 0, mesh_size*scaling_mesh_bubble};
Point(7) = {0, -radius_bubble, 0, mesh_size*scaling_mesh_bubble};
//+
Circle(1) = {4, 1, 6};
//+
Circle(2) = {6, 1, 5};
//+
Circle(3) = {5, 1, 7};
//+
Circle(4) = {7, 1, 4};
//+
Circle(5) = {2, 1, 5};
//+
Circle(6) = {5, 1, 3};
//+
Circle(7) = {3, 1, 4};
//+
Circle(8) = {4, 1, 2};
//+
Circle(9) = {2, 1, 6};
//+//+
Circle(10) = {6, 1, 3};
//+
Circle(11) = {3, 1, 7};
//+
Circle(12) = {7, 1, 2};

Point(8) = {size_rev*0.5, size_rev*0.5, size_rev*0.5, mesh_size};
Point(9) = {-size_rev*0.5, size_rev*0.5, size_rev*0.5, mesh_size};
Point(10) = {-size_rev*0.5, -size_rev*0.5, -size_rev*0.5, mesh_size};
Point(11) = {size_rev*0.5, -size_rev*0.5, -size_rev*0.5, mesh_size};
Point(12) = {size_rev*0.5, size_rev*0.5, -size_rev*0.5, mesh_size};
Point(13) = {-size_rev*0.5, -size_rev*0.5, size_rev*0.5, mesh_size};
Point(14) = {-size_rev*0.5, size_rev*0.5, -size_rev*0.5, mesh_size};
Point(15) = {size_rev*0.5, -size_rev*0.5, size_rev*0.5, mesh_size};//+
Line(13) = {9, 13};
//+
Line(14) = {13, 15};
//+
Line(15) = {15, 8};
//+
Line(16) = {8, 9};
//+
Line(17) = {9, 14};
//+
Line(18) = {14, 10};
//+
Line(19) = {10, 11};
//+
Line(20) = {11, 12};
//+
Line(21) = {12, 14};
//+
Line(22) = {8, 12};
//+
Line(23) = {11, 15};
//+
Line(24) = {10, 13};
//+
Curve Loop(1) = {10, 7, 1};
//+
//Plane Surface(1) = {-1};
Surface(1) ={-1};
//+
Curve Loop(2) = {6, -10, 2};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {6, 11, -3};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {11, 4, -7};
//+
Surface(4) = {-4};
//+
Curve Loop(5) = {9, -1, 8};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {8, -12, 4};
//+
Surface(6) = {-6};
//+
Curve Loop(7) = {12, 5, 3};
//+
Surface(7) = {-7};
//+
Curve Loop(8) = {9, 2, -5};
//+
Surface(8) = {-8};
//+
Curve Loop(9) = {22, -20, 23, 15};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {16, 17, -21, -22};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {18, 19, 20, 21};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {23, -14, -24, 19};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {14, 15, 16, 13};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {24, -13, 17, 18};
//+
Plane Surface(14) = {14};
//+
Surface Loop(1) = {13, 12, 9, 10, 14, 11};
//+
Surface Loop(2) = {1, 2, 3, 4, 6, 5, 8, 7};
//+
Volume(1) = {1, 2};
//+


Periodic Surface{11}={13} Translate {0, 0, -size_rev};
Periodic Surface{12}={10} Translate {0,-size_rev, 0};
Periodic Surface{14}={9} Translate {-size_rev,0,0};
Physical Surface(2) = {5, 6, 8, 3, 7, 2, 1, 4};
//+
Physical Surface(3) = {13, 12, 9, 10, 14, 11};
Physical Volume(1) = {1};


Mesh.ElementOrder =2;
Mesh 2;
Mesh 3;
Mesh.MshFileVersion = 2.2;