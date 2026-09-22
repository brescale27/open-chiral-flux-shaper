// Macchina Elettromeccanica 3D - Cilindro alluminio e rotore a 2 bobine
// Unità: metri
SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthMin = 0.005;
Mesh.CharacteristicLengthMax = 0.035;

r_ext = 0.05;
r_int = 0.047;
h_alu = 0.1;
r_far = 0.2;

r_coil = 0.035;
h_coil = 0.08;
r_wire = 0.006;

// 1. Rotore interno
Cylinder(1) = { 0, 0, -h_alu/2, 0, 0, h_alu, r_int };

// 2. Cilindro alluminio
Cylinder(2) = { 0, 0, -h_alu/2, 0, 0, h_alu, r_ext };
Cylinder(3) = { 0, 0, -h_alu/2, 0, 0, h_alu, r_int };
BooleanDifference(4) = { Volume{2}; Delete; }{ Volume{3}; Delete; };

// 3. Sfera Far-Field
Sphere(5) = { 0, 0, 0, r_far };

// Punti di campionamento
Point(100) = { 0.000, 0.000,  0.070 }; // Z_sopra
Point(101) = { 0.060, 0.000,  0.000 }; // XY_equatore
Point(102) = { 0.000, 0.000, -0.070 }; // Z_sotto

// Frammentazione conforme
v[] = BooleanFragments{ Volume{5}; Delete; }{ Volume{1, 4}; Point{100, 101, 102}; Delete; };

Physical Volume("AirExterior", 1) = {3};
Physical Volume("Aluminum", 2) = {2};
Physical Volume("RotorInterior", 3) = {1};
Physical Surface("FarField", 1) = {1};
