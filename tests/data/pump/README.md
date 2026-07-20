# Legacy pump regression data

`legacy_one_dimensional_reference.npz` was generated from commit `23ac736`
with the former compiled one-dimensional pump, the Tet4
`example/data/ptTet4.vtk` geometry, ASE disabled, and three 20 us
frozen-Phi-ASE RK4 steps. The archive contains every-step `betaCells`,
`betaVolume`, and `dndtPump` arrays plus JSON metadata and the geometry hash.

The general-pump regression configures one collimated 940 nm super-Gaussian
source on `ase_bottom`, a unit-transmission planar retro-relay on `ase_top`,
and total power obtained by integrating the former 16 kW/cm2 peak profile.
