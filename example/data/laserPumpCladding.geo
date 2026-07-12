SetFactory("OpenCASCADE");

// The gain volume is the 3 cm radius, 0.7 cm long laser-pump cylinder.
Cylinder(1) = {0, 0, 0, 0, 0, 0.7, 3};

// Preserve the optical boundary semantics used by the legacy slab example.
Physical Volume("gain_medium", 1) = {1};
Physical Surface("cladding", 3) = {1};
Physical Surface("ase_top", 2) = {2};
Physical Surface("ase_bottom", 1) = {3};

Mesh.CharacteristicLengthMin = 0.25;
Mesh.CharacteristicLengthMax = 0.25;
