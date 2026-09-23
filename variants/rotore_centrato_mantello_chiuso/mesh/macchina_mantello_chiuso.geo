// Macchina Elettromeccanica 3D - Variante con Mantello Chiuso a Barattolo
// Unita: metri
SetFactory("OpenCASCADE");

Mesh.Algorithm = 5;
Mesh.Algorithm3D = 1;
Mesh.CharacteristicLengthMin = 0.005;
Mesh.CharacteristicLengthMax = 0.035;

r_ext = 0.05;
r_int = 0.047;
h_alu = 0.1;
r_far = 0.2;
t_core = 0.015;
t_cap = 0.003;

// 1. Volumi rotorici interni: traferro superiore, inferiore e nucleo centrale
v_rot_up = newv; Cylinder(v_rot_up) = { 0, 0, t_core/2, 0, 0, h_alu/2 - t_core/2, r_int };
v_rot_down = newv; Cylinder(v_rot_down) = { 0, 0, -h_alu/2, 0, 0, h_alu/2 - t_core/2, r_int };
v_core = newv; Cylinder(v_core) = { 0, 0, -t_core/2, 0, 0, t_core, r_int };

// 2. Mantello cilindrico chiuso (hollow can / barattolo)
v_can_out = newv; Cylinder(v_can_out) = { 0, 0, -h_alu/2 - t_cap, 0, 0, h_alu + 2*t_cap, r_ext };
v_can_in = newv; Cylinder(v_can_in) = { 0, 0, -h_alu/2, 0, 0, h_alu, r_int };
v_mantle = newv; BooleanDifference(v_mantle) = { Volume{v_can_out}; Delete; }{ Volume{v_can_in}; Delete; };

// 3. Sfera Far-Field
v_sph = newv; Sphere(v_sph) = { 0, 0, 0, r_far };

// Frammentazione conforme
BooleanFragments{ Volume{v_sph}; Delete; }{ Volume{v_rot_up, v_rot_down, v_core, v_mantle}; Delete; };

// Assegnazione gruppi fisici
Physical Volume("AirExterior", 1) = {5};
Physical Volume("Aluminum", 2) = {4};
Physical Volume("RotorAir", 3) = {1, 2};
Physical Volume("FerromagneticCore", 4) = {3};
Physical Surface("FarField", 1) = {1};
