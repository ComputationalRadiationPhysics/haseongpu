# JuliaASE reflection surface reference

This fixture is a small single-Tet4 surface-reflection regression case.  CI does
not import or execute JuliaASE; it runs HASE on the same mesh/material data and
compares against the stored reference arrays.

`single_reflective_tet4.msh` contains one Tet4 volume and one physical surface
(id `11`) used by the surface-domain optics arrays in `reference.json`.
