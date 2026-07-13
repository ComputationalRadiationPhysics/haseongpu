# JuliaASE reflection surface reference

This fixture is a small single-Tet4 surface-reflection regression case. CI does
not import or execute JuliaASE; it runs HASE on the same mesh/material data and
compares against a high-statistics JuliaASE reference. The stored
``referenceRayCount`` is intentionally higher than CI's ``rayCount``; the latter
uses a statistical tolerance because HASE and JuliaASE use independent RNGs.

``scripts/regenerate_juliaase_reflection_fixture.py`` invokes the adjacent
JuliaASE checkout through ``scripts/juliaase_reflection_fixture.jl``. The driver
maps HASE local face 0 to JuliaASE local face 4, so physical surface 11 is the
same ``z=0`` face in both implementations. It also represents HASE's explicit
0.65 surface reflectivity as a constant JuliaASE coating, retaining total
internal reflection.

`single_reflective_tet4.msh` contains one Tet4 volume and one physical surface
(id `11`) used by the surface-domain optics arrays in `reference.json`.
